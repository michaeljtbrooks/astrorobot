#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    COORDINATE PAIRS
        
        Dimensions = single components which in combination will make up points on a grid
            Latitude
            Longitude
            HourAngle
            RightAscension
            Declination
            Azimuth
            Altitude
            
        
        Coordinatepairs = specific values on two dimensions which produces a point in space
            LatLon = a place on the Earth
            AzAlt = a place on the sky (e.g. a telescope is pointing to) in Azimuth, Altitude
            HaDec = a place on the sky (e.g. a telescope is pointing to) in HourAngle, Declination
            RADec = a place on the celestial sphere
            
            
            
        Situations:
            
            An ObservingPosition has five components:
                equatorial = (hour_angle,Declination)           #Where the equatorial mount is currently pointing (fixed unless motors move)
                celestial = (RA, Declination)                   #Where on the skymap your mount is currently pointing (varies unless tracking)
                azalt = (azimuth, altitude)                     #Where an azalt mount is currently pointing (fixed unless motors move)
                location = (latitude, longitude)                #The observer's position on the Earth
                time_travel_offset = (timestamp)                #How far forwards or back from UTC now is the observer
            
                prime_subjective = Which axis we are using to calculate all others from
            
"""
from __future__ import unicode_literals
#
from decimal import Decimal
import ephem
from LatLon import lat_lon
import math
#
from dimensions import BaseDimension, SmartLat, SmartLon, HourAngle, RightAscension, Declination, ApparentDeclination, ApparentRightAscension, Azimuth, Altitude
from libraries import sidereal
import settings
from utils import d, dms_to_hms, hms_to_dms, deg_to_rad, rad_to_deg, dd_180, coord_rotate_rad, utc_now



class BasePair(object):
    """
        Abstract class representing a point on a 2D grid defined by a pair of coordinates
    """
    X = None #This is my "x" coordinate of the pair
    Y = None
    X_name = "x" #Alternative name which can be used when fetching the coordinate
    Y_name = "y" #Alternative name which can be used when fetching the coordinate
    X_abbr = "x" #Alternative name which can be used when fetching the coordinate
    Y_abbr = "y" #Alternative name which can be used when fetching the coordinate
    X_class = BaseDimension
    Y_class = BaseDimension
    xyinverted = False #If x and y are inverted on initalisation
    name = u'' #An identifier
    apparent_position = None
    #
    def __init__(self, a=None, b=None, name=None, apparent=None):
        """
            Sets up my coordinates
        """
        if not self.xyinverted:
            x = a
            y = b
        else: #This is an inverted class!
            x = b
            y = a
        #
        #Set the name if supplied:
        if name is not None:
            self.name = unicode(name)
        #
        #Check that x and y are suitable classes
        if not isinstance(x, BaseDimension): #X isn't
            try: #Inflate the X object up into the correct class
                x = self.X_class(x)
            except TypeError as e:
                raise TypeError(u"The value for X you passed into <{my_class}> was not coercible into a <{other_class}> type coordinate".format(my_class=self.__class__.__name__, other_class=self.X_class.__name__))
        if not isinstance(y, BaseDimension): #Y isn't
            try: #Inflate the Y object up into the correct class
                y = self.Y_class(y)
            except TypeError as e:
                raise TypeError(u"The value for Y you passed into <{my_class}> was not coercible into a <{other_class}> type coordinate".format(my_class=self.__class__.__name__, other_class=self.Y_class.__name__))
        self.X = x
        self.Y = y
        
        #Deal with the "apparent" flag:
        self.apparent_position = apparent or self.apparent_position #This way mixins can override the value so long as they are put to the LEFT
    #
    def __getattr__(self, name, *args, **kwargs):
        """
        Called when there is no attribute by this name.
        We use this to grab the correct variable if one of the aliases is used
        """
        if name in (self.X_name, self.X_abbr): #Aliases for X
            return self.X
        if name in (self.Y_name, self.Y_abbr): #Aliases for X
            return self.Y
        return object.__getattribute__(self, name, *args, **kwargs)
    #
    def __setattr__(self, name, value, *args, **kwargs):
        """
        Hijacked to ensure if we set a coordinate by its alias, it ends up being stored in the right var!
        """
        if name in (self.X_name, self.X_abbr): #Aliases for X
            self.X = value
        if name in (self.Y_name, self.Y_abbr): #Aliases for X
            self.Y = value
        return object.__setattr__(self, name, value, *args, **kwargs)
    #
    def __unicode__(self):
        """
        Returns a sensible expression for these
        """
        return u"[{x_val}, {y_val}]".format(x_val=unicode(self.X), y_val=unicode(self.Y))
    #
    def __cmp__(self, other):
        """
        Compares two basepair instances. Rather senseless! Esp for wrap-around dimensions.
        
        @return: Boolean
        """
        return cmp(self.X, other.X) and cmp(self.Y, other.Y) 
    #
    def __neg__(self):
        """
        Inverse of this BasePair
        @return: BasePair instance with negative versions of X and Y (e.g. [2,3] > [-2,-3])
        """
        if self.xyinverted:
            return self.__class__(-self.Y, -self.X) #Helps us avoid flipping axes
        else:
            return self.__class__(-self.X, -self.Y) #Will work because self.X and self.Y have their own Neg function
    #
    def __pos__(self):
        if self.xyinverted:
            return self.__class__(self.Y, self.X)
        else:
            return self.__class__(self.X, self.Y)
    #
    def __abs__(self):
        self.__pos__()
        return self
    #
    def __add__(self, other, return_class=True, *args, **kwargs):
        """
        Adds two things together. Depends on what has been passed in
        
            If other is another Coordinate pair:  Extract X and Y from other, if compatible, add them
            If other is another Dimension:  If it is compatible with X or Y, add them
            If other is a tuple:  Assume t[0] is a value for self.X, and t[1] is a value for self.Y
            
            If a tuple of values supplied:
                @return: Another BasePair instance of alike class
            If another basepair class supplied:
                @return: A tuple of scalar values
        """
        if isinstance(other, self.__class__):
            #Same thing as this
            other_X = other.X
            other_Y = other.Y
            return_class = False
        elif isinstance(other, BasePair):
            #This will start to make sense for derivative classes. If something is an Instance of BaseDimension but not of the same class,
            #Then you've passed in a another type of dimension and should convert it first!
            raise TypeError(u"{my_class} object calculation expected a {my_class}, but you gave a {other_class}. Please convert it to the same type of Coordinate Pair first.".format(
                            my_class=self.__class__.__name__, other_class=other.__class__.__name__))
        elif isinstance(other, (tuple,list)):
            #Assumes (X,Y) passed in
            other_X = other[0]
            other_Y = other[1]
            if isinstance(other_X, self.X.__class__) and isinstance(other_Y, self.Y.__class__):
                #We're looking to return scalars here then!
                other_X = other_X.decimal_degree
                other_Y = other_Y.decimal_degree
                return_class = False
        elif isinstance(other, self.X.__class__):
            other_X = other.decimal_degree
            other_Y = 0
            return_class = False
        elif isinstance(other, self.Y.__class__):
            other_X = 0
            other_Y = other.decimal_degree
            return_class = False
        else:
            raise TypeError(u"{my_class} object requires an instance of {my_class}, or a coordinate pair object, or a tuple of dimension values, or a single dimension instance in order to do this calculation. You supplied a {other_class}.".format(
                            my_class=self.__class__.__name__, other_class=other.__class__.__name__))
        #
        if return_class: #We are returning a class from this
            new_X = other_X + self.X
            new_Y = other_Y + self.Y
            if self.xyinverted:
                return self.__class__(new_Y, new_X)
            else:
                return self.__class__(new_X, new_Y)
        else: #We are returning a pair of scalar values from this
            new_X = other_X + self.X.decimal_degree
            new_Y = other_Y + self.Y.decimal_degree
            return (new_X, new_Y)            
    #
    def __iadd__(self, other):
        # other is a scalar
        return self.__add__(other)
    #
    def __radd__(self, other):
        # other is a scalar
        return self.__add__(other)
    #
    def __sub__(self, other):
        """
        Subtracts something from this BasePair.
            Other coordinate pairs, dimensions, and scalars support negging.
            Tuple does not, thus we need to cater for this
        
            If a tuple or value supplied:
                @return: Another BasePair instance of alike class
            If another basepair class supplied:
                @return: A tuple of scalar values
        """
        if isinstance(other, (tuple,list)):
            neg_other = (-other[0], -other[1])
        else:
            neg_other = -other            
        return self.__add__(-neg_other)
    #
    def __isub__(self, other):
        # other is a scalar
        return self.__sub__(other)
    #
    def __rsub__(self, other):
        # other is a scalar
        return self.__sub__(other)
    #
    def __floor__(self):
        new_X = math.floor(self.X)
        new_Y = math.floor(self.Y)
        if self.xyinverted:
            return self.__class__(new_Y, new_X)
        else:
            return self.__class__(new_X, new_Y)
    #
    def __round__(self):
        new_X = round(self.X)
        new_Y = round(self.Y)
        if self.xyinverted:
            return self.__class__(new_Y, new_X)
        else:
            return self.__class__(new_X, new_Y)
    #   
    def __ceil__(self):
        new_X = math.ceil(self.X)
        new_Y = math.ceil(self.Y)
        if self.xyinverted:
            return self.__class__(new_Y, new_X)
        else:
            return self.__class__(new_X, new_Y)
    #
    def __repr__(self):
        return self.__str__()


class LatLon(BasePair, lat_lon.LatLon):
    """
    Represents a latitude longitude pair
        
        Leans on lat_lon library for clever great circle calculations
    """
    X_name = "latitude" #Alternative name which can be used when fetching the coordinate
    Y_name = "longitude" #Alternative name which can be used when fetching the coordinate
    X_abbr = "lat" #Alternative name which can be used when fetching the coordinate
    Y_abbr = "lon" #Alternative name which can be used when fetching the coordinate
    X_class = SmartLat
    Y_class = SmartLon
    #
    #Uses BasePair's init method
class LonLat(LatLon):
    xyinverted = True
    


class HourDec(BasePair):
    """
    Represents an equatorial telescope's pointing position
    """
    X_name = "hour_angle" #Alternative name which can be used when fetching the coordinate
    Y_name = "declination" #Alternative name which can be used when fetching the coordinate
    X_abbr = "ha" #Alternative name which can be used when fetching the coordinate
    Y_abbr = "dec" #Alternative name which can be used when fetching the coordinate
    X_class = HourAngle
    Y_class = ApparentDeclination
    #
    #Uses BasePair's init method
    def to_RADec(self, longitude, timestamp=None):
        """
        Converts this point to RA at the relevant earth location, at the given UTC time
        """
        #Convert X
        ra = self.X.right_ascension(longitude=longitude, timestamp=timestamp)
        dec = self.Y
        #Return new RADec object
        return RADec(ra, dec)
    def to_right_ascension(self, longitude, timestamp=None):
        """ALIAS"""
        return self.to_RADec(longitude, timestamp)
    def to_right_ascension_declination(self, longitude, timestamp=None):
        """ALIAS"""
        return self.to_RADec(longitude, timestamp)
    def to_ra(self, longitude, timestamp=None):
        """ALIAS"""
        return self.to_RADec(longitude, timestamp)
    def to_ra_dec(self, longitude, timestamp=None):
        """ALIAS"""
        return self.to_RADec(longitude, timestamp)
    #
    def to_HADec(self, *args, **kwargs):
        """To ensure standard output from this call, returns self"""
        return self
    def to_hour_angle(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_hour_angle_declination(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_ha_dec(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_ha(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    #
    def to_AzAlt(self, latitude, *args, **kwargs):
        """
        Converts this HADec into an AzAlt object provided that latitude is given
        
        @param Latitude: A latitude in decimal degrees or as a SmartLat instance 
        @return: AzAlt object for that observer
            
         For Converting HA,Dec > Az,Alt:
                x = Declination in radians (-pi to pi)
                y = Latitude in radians (-pi to pi)
                z = HourAngle in radians (0 to 2pi), (which requires longitude if converting from RA)
                
                > xt = Altitude in radians
                > yt = Azimuth in radians
            @param latitude: Observer's latitude in RADIANS
            @return HA in Radians, Dec in Radians
        """
        if not isinstance(latitude, SmartLat):
            latitude = SmartLat(latitude)
        lat_rad = latitude.radians
        ha_rad = self.X.radians
        dec_rad = self.Y.radians
        alt_rad, az_rad = coord_rotate_rad(dec_rad, lat_rad, ha_rad) #Convert to az + alt (coord rotate)
        #Cast back to Deg for conversion to an AzAlt object
        az = rad_to_deg(az_rad)
        alt = rad_to_deg(alt_rad)
        return AzAlt(az,alt)
    def to_az_alt(self, *args, **kwargs):
        return self.to_AzAlt(*args, **kwargs)
    def to_azimuth_altitude(self, *args, **kwargs):
        return self.to_AzAlt(*args, **kwargs)
    def to_az(self, *args, **kwargs):
        return self.to_AzAlt(*args, **kwargs)
class HADec(BasePair):
    pass #Alias
class DecHour(HADec):
    xyinverted = True
class DecHA(DecHour):
    pass


class RADec(BasePair):
    """
    Represents a point on the celestial sphere
    Equivalent to a specific hour angle at a given moment for a given longitude on the Earth!
    """
    X_name = "right_ascension"
    Y_name = "declination"
    X_abbr = "ra"
    Y_abbr = "dec"
    X_class = RightAscension
    Y_class = Declination
    uncorrected_position = None #Where we store our original ra_dec if using apparent
    #
    def apparent(self, observer=None, latitude=None, longitude=None, timestamp=None, temperature=settings.DEFAULT_TEMPERATURE, pressure=settings.DEFAULT_PRESSURE, elevation=settings.DEFAULT_ELEVATION):
        """
        Takes this RADec co-ordinate and maps it to the apparent position on the sky given local atmospheric refraction
        
            Uses functions from PyEphem [http://rhodesmill.org/pyephem/index.html]
        """
        if timestamp is None:
            timestamp = utc_now()
        if not observer:
            #First set up an observer
            observer = ephem.Observer()
            observer.lat = unicode(latitude.dd) #Accepts as a string only
            observer.lon = unicode(longitude.dd)#Accepts as a string only
            observer.temp = temperature
            observer.pressure = pressure
        else:
            original_observer_date = observer.date
            
        #Now update the observer to the specified time:
        observer.date = timestamp.strftime("%Y/%m/%d %H:%i:%s")
        
        if self.apparent_position: #This has already been corrected, so we need to work it all backwards to adjust for the new time / observer:
            ra_hr, ra_min, ra_sec = self.uncorrected_position.X.hms
            dec_deg, dec_min, dec_sec = self.uncorrected_position.Y.hms
        else: #Not already corrected, so use raw ra dec from self
            ra_hr, ra_min, ra_sec = self.X.hms
            dec_deg, dec_min, dec_sec = self.Y.dms
        target = ephem.readdb("Target,f|S|??,{ra_hr}:{ra_min}:{ra_sec},{dec_deg}:{dec_min}:{dec_sec},0.0,{epoch}".format(
                                ra_hr = ra_hr,
                                ra_min = ra_min,
                                ra_sec = ra_sec,
                                dec_deg = dec_deg,
                                dec_min = dec_min,
                                dec_sec = dec_sec,
                                epoch = timestamp.strftime(settings.EPHEM_DB_DATE_FORMAT)
                            ))
        target.compute(observer, epoch=observer.date) #Works out its Celestial position taking into account precession and refraction etc
        #Now generate a new object
        corrected_ra_dec = ApparentRADec(target.ra.split(":"), target.dec.split(":"))
        corrected_ra_dec.uncorrected_position = RADec(target.a_ra.split(":"), target.a_dec.split(":")) #Store the original position so we can work it back later
        #Mark it as an apparent position:
        corrected_ra_dec.apparent_position=True #Simple flag to allow us to track which positions are apparent vs real
        return corrected_ra_dec
    #
    def to_RADec(self, *args, **kwargs):
        return self
    def to_right_ascension(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    def to_right_ascension_declination(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    def to_ra(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    def to_ra_dec(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    #
    def to_HADec(self, longitude, timestamp=None):
        """
        Returns the equivalent (Hour Angle, Declination) point given a location on the earth at a particular timestamp
        
        @param longitude: a longitude as a decimal degree or SmartLon instance
        @keyword timestamp: the UTC time of the observer 
        """
        #Convert X
        ha = self.X.hour_angle(longitude=longitude, timestamp=timestamp)
        dec = self.Y
        #Create new HADec object
        return HADec(ha, dec)
    def to_hour_angle(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_hour_angle_declination(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_ha_dec(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_ha(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    #
    def to_HADec_apparent(self, observer):
        """
        Returns the equivalent (Hour Angle, Declination) for an observer on the earth with known location.
        Will correct for atmospheric refraction.
        
        @param observer: A PyEphem observer object 
        """
        ra_dec = self.apparent(observer=observer) #Whether or not self is apparent doesn't matter. The time may have changed, so need to recompute it
        #Convert X
        ha = ra_dec.X.hour_angle(longitude=observer.lon, observer=observer.date)
        dec = ra_dec.Y
        #Create new HADec object
        return HADec(ha, dec)
    def to_hour_angle_apparent(self, *args, **kwargs):
        return self.to_HADec_apparent(*args, **kwargs)
    def to_hour_angle_declination_apparent(self, *args, **kwargs):
        return self.to_HADec_apparent(*args, **kwargs)
    def to_ha_dec_apparent(self, *args, **kwargs):
        return self.to_HADec_apparent(*args, **kwargs)
    def to_ha_apparent(self, *args, **kwargs):
        return self.to_HADec_apparent(*args, **kwargs)
    #
    def to_AzAlt(self, latitude, longitude, timestamp=None, *args, **kwargs):
        """
        Converts this RADec into an AzAlt object provided that latitude, longitude and timestamp is given
        
        @param Latitude: A latitude in decimal degrees or as a SmartLat instance
        @param Longitude: A longitude in decimal degrees or as a SmartLon instance 
        @return: AzAlt object for that observer
        """
        #First convert to Hour Angle, Dec object:
        ha_dec = self.to_HADec(longitude, timestamp) #Will be corrected for refraction here!
        #Then convert the HADec object to AzAlt
        az_alt = ha_dec.to_AzAlt(latitude)
        return az_alt
    def to_az_alt(self, *args, **kwargs):
        return self.to_AzAlt(*args, **kwargs)
    def to_azimuth_altitude(self, *args, **kwargs):
        return self.to_AzAlt(*args, **kwargs)
    def to_az(self, *args, **kwargs):
        return self.to_AzAlt(*args, **kwargs)
class RightascensionDeclination(RADec):
    pass #Alias
class RightAscensionDeclination(RADec):
    pass #Alias
class DecRA(RADec):
    xyinverted = True
class DeclinationRightAscension(DecRA):
    pass
class DeclinationRightascension(DecRA):
    pass

class ApparentRADec(RADec):
    """
    A variant of RADec showing an "apparent" position
    """
    apparent_position=True
#ALIASES
ApparentRightascensionDeclination = ApparentRADec
ApparentRightAscensionDeclination = ApparentRADec
class ApparentDecRA(ApparentRADec):
    xyinverted = True
ApparentDeclinationRightAscension = ApparentDecRA
ApparentDeclinationRightascension = ApparentDecRA


class AzAlt(BasePair):
    """
    Represents the location relative to the local geography that a scope is pointing to
        
        Given a known latitude, can map to HourAngle,Dec
        Given a known latitude+longitude+time, can map to RA,Dec
    """
    X_name = "azimuth"
    Y_name = "altitude"
    X_abbr = "az"
    Y_abbr = "alt"
    X_class = Azimuth
    Y_class = Altitude
    #
    def to_RADec(self, latitude, longitude, timestamp=None):
        """
        Maps to (Right Ascension, Declination)
        
            First casts to HA,Dec using latitude, 
            Then RA,Dec using longitude+timestamp
        """
        #First, get HADec object:
        ha_dec_obj = self.to_ha_dec(latitude)
        #Now convert the HADec instance into RADec
        ra_dec_obj = ha_dec_obj.to_ra_dec(longitude, timestamp)
        return ra_dec_obj
    def to_right_ascension(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    def to_right_ascension_declination(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    def to_ra(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    def to_ra_dec(self, *args, **kwargs):
        """ALIAS"""
        return self.to_RADec(*args, **kwargs)
    #
    def _rotate_to_HADec(self, lat_rad):
        """
        Returns the equivalent (Hour Angle, Declination) point for your AzAlt at your latitude
            
            For Converting Az,Alt > RA,Dec:
                x = Altitude in radians (-pi to pi)
                y = Latitude in radians (-pi to pi)
                z = Azimuth in radians (0 to 2pi)
                
                > xt = Declination in radians
                > yt = HourAngle in radians
        
            @param lat_rad: Observer's latitude in RADIANS
            @return HA in Radians, Dec in Radians
        """
        az_rad = self.X.radians
        alt_rad = self.Y.radians
        dec_rad, ha_rad = coord_rotate_rad(alt_rad, lat_rad, az_rad) #Convert to hour angle + dec (coord rotate)
        #Now return in correct order
        return (ha_rad, dec_rad,)
    #
    def to_HADec(self, latitude, *args, **kwargs):
        """
        Returns the equivalent (Hour Angle, Declination) point for your AzAlt at your latitude
            
        @param latitude: The observer's latitude in decimal degrees or as a SmartLat object
        @return: A HADec coordinate pair object 
        """
        #Ensure latitude is in correct format
        if not isinstance(latitude, SmartLat): #Cast strings / decimals for longitude into a SmartLon Longitude object
            latitude = SmartLat(Decimal(latitude))
        #Convert to radians
        lat_rad = latitude.radians
        ha_rad, dec_rad = self._rotate_to_HADec(lat_rad)
        ha = rad_to_deg(ha_rad)
        dec = rad_to_deg(dec_rad)
        #Return new HADec object
        return HADec(ha, dec)
    def to_hour_angle(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_hour_angle_declination(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_ha_dec(self, *args, **kwargs):
        return self.to_HADec(*args, **kwargs)
    def to_ha(self, *args, **kwargs): #Technically a misnomer!
        return self.to_HADec(*args, **kwargs)
    #
    def to_AzAlt(self, *args, **kwargs):
        return self #Already is!
    def to_az_alt(self, *args, **kwargs):
        return self #Already is!
    def to_azimuth_altitude(self, *args, **kwargs):
        return self #Already is!
    def to_az(self, *args, **kwargs):
        return self
#
#ALIASES
class AzimuthAltitude(AzAlt):
    pass
class AltAz(AzAlt):
    xyinverted = True
class AltitudeAzimuth(AltAz):
    pass

