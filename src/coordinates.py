#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    Coordinates
        
        This is the representation of where objects are on the Earth and in space. Organised into three tiers of base objects:
        
            coordinate = one component of a pair of coordinates. Has all the machinery to convert between Deg / Rad / Hrs
            gridlocation = a pair of coordinates
        
        
        Some standards to adopt:
        
            A SkyPosition has three components:
                equatorial = (hour_angle,Declination)           #Where to point your equatorial mount (varies with time)
                    - hour_angle is in HOURS
                celestial = (RA, Declination)                   #Where object is on the skymap (fixed)
                    - RA is in HOURS
                azalt = (azimuth, altitude)                     #Where to point your AzAlt mount (varies with time)
                    - azimuth is in DEGREES
            
            A MountPosition has three components too:
                equatorial = (hour_angle,Declination)           #Where the equatorial mount is currently pointing (fixed unless motors move)
                celestial = (RA, Declination)                   #Where on the skymap your mount is currently pointing (varies unless tracking)
                azalt = (azimuth, altitude)                     #Where an azalt mount is currently pointing (fixed unless motors move)
                
            Thus, current time and current location are necessary params for working out 
            
    
    @TODO: Banish libraries which are using floats! Rewrite to use decimals
        
"""

from datetime import datetime, timedelta
from decimal import Decimal 
import math
#
from dateutil import parser
from LatLon import lat_lon, LatLon, Latitude, Longitude
from libraries import sidereal
import pytz


### Constants ###
PI = Decimal("3.141592653589793238462643383279502884197169399") #I hate floats! This is a far more accurate Decimal Pi
TWO_PI = 2 * PI


def d(val, *args, **kwargs):
    """
    Casts the value into a true decimal
    """
    val_str = unicode(val)
    val_d = Decimal(val_str, *args, **kwargs)
    return val_d


def hms_to_dms(hours=Decimal("0.0"), minutes=Decimal("0.0"), seconds=Decimal("0.0")):
    """
    Converts a time based coordinate to a degrees based one
    """
    #Cast to Decimal for accuracy:
    seconds = Decimal(unicode(seconds))
    minutes = Decimal(unicode(minutes))
    hours = Decimal(unicode(hours)) 
    #
    #Degrees are 15x hours:
    deg_seconds_decimal = seconds * Decimal("15.0")
    deg_minutes_carryover, deg_seconds = divmod(deg_seconds_decimal, Decimal("60"))
    deg_minutes_decimal = minutes * Decimal("15.0") + deg_minutes_carryover
    deg_hours_carryover, deg_minutes = divmod(deg_minutes_decimal, Decimal("60"))
    deg_hours = hours * Decimal("15.0") + deg_hours_carryover
    #
    return(deg_hours,deg_minutes,deg_seconds)


def dms_to_hms(degrees=Decimal("0.0"), minutes=Decimal("0.0"), seconds=Decimal("0.0")):
    """
    Converts a degrees based coordinate into a time based one
    """
    #Cast to Decimal for accuracy:
    deg_seconds = Decimal(unicode(seconds))
    deg_minutes = Decimal(unicode(minutes))
    deg_degrees = Decimal(unicode(degrees))
    #
    #Hours are Degrees/15
    hr_hours, deg_degrees_remainder = divmod(deg_degrees, Decimal("15.0"))
    #Minutes are DegMinutes/15 
    deg_minutes_carrydown = deg_degrees_remainder * Decimal("60.0")
    hr_minutes, deg_minutes_remainder = divmod((deg_minutes + deg_minutes_carrydown), Decimal("15.0"))
    #Seconds are DegSeconds/15
    deg_seconds_carrydown = deg_minutes_remainder * Decimal("60.0")
    hr_seconds = ((deg_seconds + deg_seconds_carrydown)/Decimal("15.0")) #Last entity can be a decimal
    #
    return(hr_hours, hr_minutes, hr_seconds)
    
def dd_180(degrees=Decimal("0")):
    """
    Converts decimal degrees in range 180-360 to a negative equivalent
    """
    degrees = Decimal(unicode(degrees)) #Cast to decimal
    out_degrees = degrees % Decimal("360")
    if out_degrees > Decimal("180"):
        out_degrees = out_degrees - Decimal("360")  
    return out_degrees

def utc_now():
    """
    Returns current true UTC time as a time-aware datetime object
    """
    realtime = datetime.utcnow()
    realtime = pytz.utc.localize(realtime)
    return realtime

def deg_to_rad(degrees):
    """
    Purist Decimal method to convert degrees to radians, normalised to 0-2pi
    """
    degrees = Decimal(unicode(degrees))
    radians = degrees * TWO_PI / 360
    radians = radians % (TWO_PI)
    return radians

def rad_to_deg(radians):
    """
    Purist Decimal method for radians to degrees
    """
    radians = Decimal(unicode(radians))
    degrees = radians * (360/TWO_PI)
    degrees = degrees % (Decimal("360"))
    return degrees


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
    #
    def __init__(self, big=Decimal("0"), medium=Decimal("0"), small=Decimal("0"), mode=None):
        """
            Sets up this object based upon its format.
            
            @TODO: handle degrees in range 180-360 as negative degrees (-180-0)
        """
        #Defaults
        if mode:
            if self.mutable_mode:
                self.base_mode = mode
        else:
            mode = self.base_mode
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
            return u"{0}h {1}m {2:.5f}s".format(hours,minutes,seconds)
        elif mode == "hd": #Hours as decimal:
            return u"{0}h".format(self.hd)
        elif mode == "rad": #Radians
            return u"{0} radians".format(self.rad)
        elif mode == "dms": #Degrees min secs
            return u"{0}° {1}' {2:.5f}\"".format(*self.dms)
        elif mode == "dd": #Decimal degrees
            return u"{0}°".format(self.dd)
        return unicode(self.decimal_degree)
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
    """ 
    base_mode = "hms"
    mutable_mode = False
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


class Declination(BaseCoordinate):
    """
    ##HERE##
    """










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
        
    
