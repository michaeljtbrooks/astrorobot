#!/usr/bin/python
# -*- coding: utf8 -*-

"""


OBSOLETE - use dimensions.py instead!!!

OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE OBSOLETE 
        
"""

from datetime import datetime, timedelta
from dateutil import parser
from decimal import Decimal 
from LatLon import lat_lon
import math
import pytz
#
from libraries import sidereal
from utils import d, dms_to_hms, hms_to_dms, deg_to_rad, rad_to_deg, dd_180, coord_rotate_rad, utc_now


################### Classes #########################

#### Single points

class BaseCoordinate(lat_lon.GeoCoord, object):
    """
    Provides methods for converting between:
    
        dd     Decimal degrees
        dms    Degrees, minutes, seconds
        hms    Hours, minutes, seconds
        hd     Decimal hours
        rad    Radian

        Depends upon LatLon's GeoCoord methods
    """        
    
    base_mode = "dms" #Whether to assume this coordinate operates in degrees / hours or radians
    mutable_mode = True #Whether we can piss about with the mode and change it once initialised
    str_prefix = None
    range_high = 360
    range_low = 0
    range_rotates = True #True = The range flicks from high back to low if you increase. False = Range oscillates  
    #
    def __init__(self, big=Decimal("0"), medium=Decimal("0"), small=Decimal("0"), mode=None):
        """
            Sets up this object based upon its format.
            
            @TODO: handle a range of string inputs!
            @TODO: handle degrees in range 180-360 as negative degrees (-180-0)
        """
        #Defaults
        if mode:
            if self.mutable_mode:
                self.base_mode = mode
        else:
            mode = self.base_mode
        #Parsing input
        if isinstance(big,(tuple,list)):
            val_tup = big
            try:
                big = val_tup[0]
            except IndexError: #Big number mandatory!
                raise TypeError("You passed in an empty list or tuple to %s. Need at least one value!" % self.__class__.__name__)
            try:
                medium = val_tup[1]
                small = val_tup[2]
            except IndexError:
                pass #These smaller entities are not mandatory!
        #Setting to decimal
        big = Decimal(unicode(big))
        medium = Decimal(unicode(medium))
        small = Decimal(unicode(small))
        degrees = big
        minutes = medium
        seconds = small
        #Handling formats
        if mode == "hms": #Hours, minutes, seconds supplied
            self.hr_hours = big
            self.hr_minutes = medium
            self.hr_seconds = small
            degrees, minutes, seconds = hms_to_dms(big, medium, small)
        elif mode == "hd": #Hours as decimal degrees:
            degrees = big * Decimal("15.0")
        elif mode == "rad": #Radians
            degrees = Decimal(math.degrees(big)) % Decimal("360")
        #Now we can pass this into super
        super(BaseCoordinate, self).__init__(degrees, minutes, seconds)
    #
    def normalised(self, val=None, as_float=False):
        """
        Returns the specified value to fit the prescribed ranges and range rules
        """
        if val is None:
            val = self.decimal_degree
        val = d(val) #Cast to Decimal
        range_high = d(self.range_high)
        range_low = d(self.range_low)
        if self.range_rotates:
            #Cycles through the range
            n_val = (val + range_high)%(range_high - range_low) + range_low
        else:
            #Oscillates between range extremes:
            range_range = range_high - range_low
            n_passes, residual = divmod((val - range_low),range_range)
            if n_passes % 2: #If odd, subtract residual off the high
                n_val = range_high - residual
            else: #If even, add it to the low
                n_val = range_low + residual
        return n_val                 
    #
    @property
    def hd(self):
        """Returns value in decimal hours"""
        return Decimal(self.decimal_degree)/Decimal("15.0")
    #
    @property
    def hms(self):
        """Returns value in hours, minutes and seconds - useful for RA & hour_angle calculations"""
        self.hr_hours, self.hr_minutes, self.hr_seconds = dms_to_hms(self.degree, self.minute, self.second)
        return(self.hr_hours, self.hr_minutes, self.hr_seconds)
    #
    @property
    def dd(self):
        """Returns value in decimal degrees"""
        return self.decimal_degree
    #
    @property
    def dd360(self):
        """Returns value in decimal degrees, normalised to range 0-360"""
        return self.decimal_degree % Decimal("360")
    @property
    def range360(self):
        return self.dd360
    #
    @property
    def dd180(self):
        """Returns the decimal value adjusted to range -180 - 180"""
        out_dd180 = dd_180(self.decimal_degree)
        return out_dd180
    @property
    def range180(self):
        return self.dd180    
    #
    @property
    def dms(self):
        """Returns value as tuple of degrees, minutes, seconds"""
        return (self.degree, self.minute, self.second)
    #
    @property
    def hemisphere_dd(self):
        """
        A string expression of the coordinate with hemisphere applied 
        """
        return "{dd}{hemisphere}".format(dd=self.decimal_degree, hemisphere=self.get_hemisphere())
    #
    @property
    def hemisphere_dms(self):
        """
        A string expression of the coordinate as deg min sec with hemisphere applied (if one exists) 
        """
        if self.get_hemisphere():
            out = u"{0} {1} {2:.5f}{hemisphere}".format(self.degree, self.minute, self.second, hemisphere=self.get_hemisphere())
        else:
            out = self._signed_unicode(u"{0} {1} {2:.5f}{hemisphere}", self.degree, self.minute, self.second, hemisphere=self.get_hemisphere())
        if self.str_prefix:
            out = u"{0} {1}".format(unicode(self.str_prefix), out)
        return out
    #
    @property
    def radians(self):
        """
        Gives the co-ordinates' position in radians for use with sidereal.py calculations
        
        @return: This coordinate's position in radians in range 0 to 2π
        """
        coord_str = self.hemisphere_dd
        if self.get_hemisphere in ("N","S"):
            angle = sidereal.parseLat(coord_str)
        elif self.get_hemisphere in ("E","W"):
            angle = sidereal.parseLon(coord_str)
        else:
            angle = self._calc_radians()
        return angle
    @property
    def rad(self):
        """ALIAS"""
        return self.radians
    #
    def _calc_radians(self):
        """
        Calculates the radians in range 0 to 2π from given decimal degrees
        """
        radian = deg_to_rad(self.decimal_degree)
        return radian
    #
    def get_hemisphere(self):
        """
        Override this as parent class does mandates it to be derived as it is an abstract class
        """
        return ''
    #
    def set_hemisphere(self, format_str):
        """
        Override this as parent class does mandates it to be derived as it is an abstract class
        """
        return ''
    #
    def _signed_unicode(self, format_str, *args, **kwargs):
        """
        Formats a tuple as a string with a prefixed sign
        
        @param format_str: A string or unicode defining the output format
        @args: the data that will be flushed into the format string 
        """
        sign = u"+"
        if args[0] < 0:
            sign = u"-"
        #Now flip all to positive
        new_args = []
        for arg in args:
            new_arg = arg
            if arg < 0:
                new_arg = 0 - arg
            new_args.append(new_arg)
        out_unicode = sign + unicode(format_str.format(*new_args, **kwargs))
        return out_unicode        
    #
    def __unicode__(self):
        """
            Returns a string representation of this, depending on its format
        """
        mode = self.base_mode
        if mode == "hms": #Hours, minutes, seconds supplied
            hours, minutes, seconds = self.hms
            out = u"{0}h {1}m {2:.5f}s".format(hours,minutes,seconds)
        elif mode == "hd": #Hours as decimal:
            out = u"{0}h".format(self.hd)
        elif mode == "rad": #Radians
            out = u"{0} radians".format(self.rad)
        elif mode == "dms": #Degrees min secs
            out = u"{0}° {1}' {2:.5f}\"".format(*self.dms)
        elif mode == "dd": #Decimal degrees
            out = u"{0}°".format(self.dd)
        else:
            out = unicode(self.decimal_degree)
        if self.str_prefix:
            out = u"{0} {1}".format(unicode(self.str_prefix), out)
        return out
    #
    def __str__(self):
        """
        Returns what unicode returns, but sans any special chars
        """
        out = self.__unicode__()
        out = out.replace(u"°",u"*")
        return str(out)
    #
    def type(self):
        """
        Identifies the object type
        """
        return self.__class__.__name__


class SmartLat(BaseCoordinate):
    """
    Represents an earth Latitude, with full casting capabilities
        
        Derived from: GeoCoord by Gen Del Raye
    """
    base_mode = "dms"
    mutable_mode = False
    str_prefix = "Lat"
    range_high = 90
    range_low = -90
    range_rotates = False #True = The range flicks from high back to low if you increase. False = Range oscillates
    
    def range90(self):
        """
        Reports the latitude in range -90 to +90
        """
    
    def get_hemisphere(self):
        '''
        Returns the hemisphere identifier for the current coordinate
        '''
        if self.decimal_degree < 0: return 'S'
        else: return 'N'
        
    def set_hemisphere(self, hemi_str):
        '''
        Given a hemisphere identifier, set the sign of the coordinate to match that hemisphere
        '''
        if hemi_str == 'S':
            self.degree = abs(self.degree)*-1
            self.minute = abs(self.minute)*-1
            self.second = abs(self.second)*-1
            self._update()
        elif hemi_str == 'N':
            self.degree = abs(self.degree)
            self.minute = abs(self.minute)
            self.second = abs(self.second)
            self._update()
        else:
            raise ValueError('Hemisphere identifier for latitudes must be N or S')
        
    def __repr__(self):
        return 'Latitude %s' %(self.__str__())
class Latitude(SmartLat):
    """ALIAS"""
    pass


class SmartLon(BaseCoordinate):
    """
    Represents an earth Longitude, with full casting capabilities
    
        Derived from: GeoCoord by Gen Del Raye
    """
    base_mode = "dms"
    mutable_mode = False
    str_prefix = "Lon"

    def __init__(self, degree = 0, minute = 0, second = 0):
        super(SmartLon, self).__init__(degree, minute, second) # Initialize the GeoCoord
        decimal_degree = self.range180() # Make sure that longitudes are reported in the range -180 to 180
        self.degree, self.minute, self.decimal_minute, self.second = self._calc_degreeminutes(decimal_degree)
        self._update()
    
    def range180(self):
        '''
        Report longitudes using the range -180 to 180.
        '''
        return ((self.decimal_degree + 180)%360) - 180
        
    def range360(self):
        '''
        Report longitudes using the range 0 to 360
        '''
        return (self.decimal_degree + 360)%360
        
    def get_hemisphere(self):
        '''
        Returns the hemisphere identifier for the current coordinate
        '''
        if self.decimal_degree < 0: return 'W'
        else: return 'E'

    def set_hemisphere(self, hemi_str):
        '''
        Given a hemisphere identifier, set the sign of the coordinate to match that hemisphere
        '''
        if hemi_str == 'W':
            self.degree = abs(self.degree)*-1
            self.minute = abs(self.minute)*-1
            self.second = abs(self.second)*-1
            self._update()
        elif hemi_str == 'E':
            self.degree = abs(self.degree)
            self.minute = abs(self.minute)
            self.second = abs(self.second)
            self._update()
        else:
            raise(ValueError, 'Hemisphere identifier for longitudes must be E or W')
    
    def __repr__(self):
        return 'Longitude %s' %(self.__str__())
class Longitude(SmartLon):
    """ALIAS"""
    pass


class HourAngle(BaseCoordinate):
    """
    The horizontal pointing of an equatorial mounted scope.
        
        Definitions
            Local Hour Angle = 0 = is where the scope is pointing to the plane which passes through the Earth's axis and the Zeneth
            Greenwich Hour Angle = 0 = is where the scope is pointing to the plane which passes through Longitude 0
        
        When converting between RA and LA, this uses LOCAL HOUR ANGLE
    """
    base_mode = "hms"
    mutable_mode = False
    str_prefix = "HrA"
    range_high = 24
    range_low = 0
    range_rotates = True #True = The range flicks from high back to low if you increase. False = Range oscillates
    #
    def right_ascension(self, longitude, timestamp=None):
        """
        Gives the Right Ascension of this hour_angle coordinate if time and longitude of observer supplied!
        
        @param longitude: The longitude of the observer (as a true SmartLon object, or a string)
        @keyword timestamp:    The UTC time of the observer (will assume now if none supplied)
        
        @returns an RightAscension object with the RA in it!
        """
        if timestamp is None: #Default to real now
            timestamp = utc_now()
        if not isinstance(longitude, SmartLon): #Cast strings / decimals for longitude into a SmartLon Longitude object
            longitude = SmartLon(Decimal(longitude))
        #Convert all our values into radian floats which sidereal loves:
        long_rad = float(longitude.radians) 
        hr_rad = float(self.radians)
        #Now call upon Sidereal's function
        ra_rad =  sidereal.hourAngleToRA(hr_rad, timestamp, long_rad) #Outputs hour angle in radians
        #Now inflate the output hour angle into a proper object
        return RightAscension(ra_rad, mode="rad") #Will accept a rad, and return a properly filled HourAngle
    def to_right_ascension(self, *args, **kwargs):
        """ALIAS"""
        return self.right_ascension(*args, **kwargs)
    def to_ra(self, *args, **kwargs):
        """ALIAS"""
        return self.right_ascension(*args, **kwargs)
    #
    def hour_angle(self, *args, **kwargs):
        """
        This is already an hour_angle, so return self!
        """
        return self
    def to_hour_angle(self, *args, **kwargs):
        return self
    def to_ha(self, *args, **kwargs):
        return self


class RightAscension(BaseCoordinate):
    """
    Represents the Right Ascension location of an object on the celestial sphere
    
        0,24      Vernal Equinox
        12     Opposite Vernal Equinox
    """ 
    base_mode = "hms"
    mutable_mode = False
    str_prefix = "RA"
    range_high = 24
    range_low = 0
    range_rotates = True #True = The range flicks from high back to low if you increase. False = Range oscillates
    #
    def hour_angle(self, longitude, timestamp=None):
        """
        Gives the hour_angle of this RA coordinate if time and longitude of observer supplied!
        
        @param longitude: The longitude of the observer (as a true SmartLon object, or a string)
        @keyword timestamp:    The UTC time of the observer (will assume now if none supplied)
        
        @returns an HourAngle object with the HourAngle in it!
        """
        if timestamp is None: #Default to real now
            timestamp = utc_now()
        if not isinstance(longitude, SmartLon): #Cast strings / decimals for longitude into a SmartLon Longitude object
            longitude = SmartLon(Decimal(longitude))
        #Convert all our values into radian floats which sidereal loves:
        long_rad = float(longitude.radians) 
        ra_rad = float(self.radians)
        #Now call upon Sidereal's function
        hour_angle_rad =  sidereal.raToHourAngle(ra_rad, timestamp, long_rad) #Outputs hour angle in radians
        #Now inflate the output hour angle into a proper object
        return HourAngle(hour_angle_rad, mode="rad") #Will accept a rad, and return a properly filled HourAngle
    def to_hour_angle(self, *args, **kwargs):
        """ALIAS"""
        return self.hour_angle(*args, **kwargs)
    def to_ha(self, *args, **kwargs):
        """ALIAS"""
        return self.hour_angle(*args, **kwargs)
    #
    def right_ascension(self, *args, **kwargs):
        """
        This is already a RA, so return self!
        """
        return self
    def to_right_ascension(self, *args, **kwargs):
        return self
    def to_ra(self, *args, **kwargs):
        return self


class Declination(SmartLat):
    """
    Represents the Declination (N+ S- direction) which applies to both the celestial grid and equatorial scopes (they are aligned!!)
    
        [Also almost perfectly aligned to Latitude, but see: https://en.wikipedia.org/wiki/Declination#Relation_to_latitude
         since the two are separate concepts we will use two separate classes]
         
        +90    Celestial North pole
        0      Celestial equator
        -90    Celestial South pole
        
    """
    base_mode = "dms"
    mutable_mode = False
    str_prefix = "Dec"
    range_high = 90
    range_low = -90
    range_rotates = False #True = The range flicks from high back to low if you increase. False = Range oscillates
    

class Altitude(BaseCoordinate):
    """
    Represents the Altitude angle of an AltAz scope. +90 deg is the Zeneth, and -90 deg is the floor and utterly pointless to point to.
    Remember, the horizon is slightly below 0 deg for flat earth as the earth is curved and telescopes have some height
    
        +90    Zeneth
        0      Geometrical Horizon
        -90    Nadir
    """
    base_mode = "dms"
    mutable_mode = False
    str_prefix = "Alt"
    range_high = 90
    range_low = -90
    range_rotates = False #True = The range flicks from high back to low if you increase. False = Range oscillates
    
    #You cannot convert from Altitude alone to Declination because it could be any of a range of values depending on the Azimuth!
    

class Azimuth(BaseCoordinate):
    """
    Represents the Azimuth angle of an AltAz scope. Range 0 to 360,  
    
        0,360  North
        180    South
    """
    base_mode = "dms"
    mutable_mode = False
    str_prefix = "Alt"
    range_high = 360
    range_low = 0
    range_rotates = True #True = The range flicks from high back to low if you increase. False = Range oscillates



#### Pair of coordinates

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
    X_class = BaseCoordinate
    Y_class = BaseCoordinate
    xyinverted = False #If x and y are inverted on initalisation
    name = u'' #An identifier
    #
    def __init__(self, a=None, b=None, name=None):
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
        if not isinstance(x, BaseCoordinate): #X isn't
            try: #Inflate the X object up into the correct class
                x = self.X_class(x)
            except TypeError as e:
                raise TypeError("The value for X you passed into <{my_class}> was not coercible into a <{other_class}> type coordinate".format(my_class=self.__class__.__name__, other_class=self.X_class.__name__))
        if not isinstance(y, BaseCoordinate): #Y isn't
            try: #Inflate the Y object up into the correct class
                y = self.Y_class(y)
            except TypeError as e:
                raise TypeError("The value for Y you passed into <{my_class}> was not coercible into a <{other_class}> type coordinate".format(my_class=self.__class__.__name__, other_class=self.Y_class.__name__))
        self.X = x
        self.Y = y
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
    X_class = SmartLat
    Y_class = SmartLon
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
        """
        #Convert X
        ha = self.X.hour_angle(longitude=longitude, timestamp=timestamp)
        dec = self.Y
        #Return new RADec object
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
class RightascensionDeclination(RADec):
    pass #Alias
class RightAscensionDeclination(RADec):
    pass #Alias


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
        
            @param latitude: Observer's latitude in RADIANS
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


class EarthLocation(LatLon):
    """
    Represents a position on the Earth's surface. Allows you to return all sorts of useful things given your location:
    
        Location as Lat/Lon
        Location as Lat Radians / Lon Radians
    """
    time_travel_offset = None
    #
    def __init__(self, lat='0', lon='0', time=None, *args, **kwargs):
        """
        @param lat: Your position's Latitude as a STRING in DEGREES
        @param lon: Your position's Longitude as a STRING in DEGREES
        
        @keyword time: The time you wish to emulate. Assumes current UTC time unless you say otherwise  
        
        """
        #Deal with input coordinates
        if not (isinstance(lat, str) and isinstance(lon, str)):
            raise SyntaxError("AstroRobot > EarthLocation: Please specify your latitude and longitude as two strings e.g lat='51.456', lon='-0.721'")
        lat_obj = SmartLat(Decimal(lat))
        lon_obj = SmartLon(Decimal(lon))
        #
        #Deal with time:        
        if time is not None:
            self.timetravel(time)
        else: #No time provided, so offset is zero and default to utc_now
            self.time_travel_offset = timedelta(seconds=0)
        #
        #Inflate into full object
        super(EarthLocation, self).__init__(lat=lat_obj, lon=lon_obj, *args, **kwargs)
    #
    @property
    def utc_real(self):
        """
        Returns current true UTC time
        """
        realtime = datetime.utcnow()
        realtime = pytz.utc.localize(realtime)
        return realtime
    #
    @property
    def now(self):
        """
        Returns the current time after timetravel has been applied!
        
        @return: UTC timezone-aware "current" time, after timetravel
        """
        return self.utc_real + self.time_travel_offset
    #
    def timetravel(self, time):
        """
        Time travels you to the specified time!
        
        @param time: The time you wish to emulate.
        
        @return: A datetime object for your current time
        """
        if not isinstance(time, datetime): #Convert whatever has been supplied into a datetime
            time = parser.parse(unicode(time)) #Assume string
            time = pytz.utc.localize(time)
        #Calculate offset
        self.time_travel_offset = time - self.utc_real
        return time #Returns the correctly inflated datetime object representing now
    #
    @property
    def radians(self):
        """
        Returns the coordinates Lat/Lon in RADIANS (for use with sidereal!)
        """
        return (self.lat.radians, self.lon.radians)
    #
    def convert_ra_to_hourangle(self, ra):
        """
        Converts a given Right Ascension (in Degrees)
        """
        
    
