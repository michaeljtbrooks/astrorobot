#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    DIMENSIONS
        
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

from decimal import Decimal 
from LatLon import lat_lon
import math
#
from libraries import sidereal
from utils import d, dms_to_hms, hms_to_dms, deg_to_rad, dd_180, utc_now


class BaseDimension(object):
    """
    Provides methods for converting between:
    
        dd     Decimal degrees
        dms    Degrees, minutes, seconds
        hms    Hours, minutes, seconds
        hd     Decimal hours
        rad    Radian

        Adapted from LatLon by Gen Del Raye
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
        #Now we can set these as class objects
        self.degree = degrees
        self.minute = minutes
        self.second = seconds
        self._update() # Clean up each variable and make them consistent

    #
    def _update(self):
        '''
        Given degree, minute, and second information, clean up the variables and make them
        consistent (for example, if minutes > 60, add extra to degrees, or if degrees is
        a decimal, add extra to minutes).
        '''
        self.decimal_degree = self._calc_decimaldegree(self.degree, self.minute, self.second)
        self.degree, self.minute, self.decimal_minute, self.second = self._calc_degreeminutes(self.decimal_degree)
    #
    def set_minute(self, minute):
        self.minute = d(minute)
    #    
    def set_second(self, second):
        self.second = d(second)
    #
    def set_degree(self, degree):
        self.degree = d(degree)
    #
    @staticmethod
    def _calc_decimaldegree(degree, minute, second):
        '''
        Calculate decimal degree form degree, minute, second
        '''
        return d(degree) + d(minute)/d("60.0") + d(second)/d("3600.0")
    
    @staticmethod
    def _calc_degreeminutes(decimal_degree):
        '''
        Calculate degree, minute second from decimal degree
        '''
        sign = cmp(decimal_degree, 0) # Store whether the coordinate is negative or positive
        decimal_degree = abs(decimal_degree)
        degree = decimal_degree//1 # Truncate degree to be an integer
        decimal_minute = (decimal_degree - degree)*60. # Calculate the decimal minutes
        minute = decimal_minute//1 # Truncate minute to be an integer
        second = (decimal_minute - minute)*60. # Calculate the decimal seconds
        # Finally, re-impose the appropriate sign
        degree = degree*sign
        minute = minute*sign
        second = second*sign
        return (degree, minute, decimal_minute, second)
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
    def to_string(self, format_str):
        '''
        Output lat, lon coordinates as string in chosen format
        Inputs:
            format (str) - A string of the form A%B%C where A, B and C are identifiers.
              Unknown identifiers (e.g. ' ', ', ' or '_' will be inserted as separators 
              in a position corresponding to the position in format.
        Examples:
            >> palmyra = LatLon(5.8833, -162.0833)
            >> palmyra.to_string('D') # Degree decimal output
            ('5.8833', '-162.0833')
            >> palmyra.to_string('H% %D')
            ('N 5.8833', 'W 162.0833')
            >> palmyra.to_string('d%_%M')
            ('5_52.998', '-162_4.998')
        '''
        format2value = {'H': self.get_hemisphere(),
                        'M': abs(self.decimal_minute),
                        'm': int(abs(self.minute)),
                        'd': int(self.degree),
                        'D': self.decimal_degree,
                        'S': abs(self.second)}
        format_elements = format_str.split('%')
        coord_list = [str(format2value.get(element, element)) for element in format_elements]
        coord_str = ''.join(coord_list)
        if 'H' in format_elements: # No negative values when hemispheres are indicated
            coord_str = coord_str.replace('-', '')
        return coord_str
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
    #
    def __cmp__(self, other):
        """
        Compares two dimension instances
        """
        if isinstance(other, self.__class__):
            #Same thing as 
            other = other.decimal_degree
        elif isinstance(other, BaseDimension):
            #This will start to make sense for derivative classes. If something is an Instance of BaseDimension but not of the same class,
            #Then you've passed in a another type of dimension and should convert it first!
            raise TypeError(u"{my_class} object calculation expected a {my_class}, but you gave a {other_class}. Please convert it to the same type of dimension first.".format(
                            my_class=self.__class__.__name__, other_class=other.__class__.__name__)
                            )
        else:
            other = d(other) #Otherwise, try to cast it to a Decimal
        return cmp(self.decimal_degree, other.decimal_degree)
    #
    def __neg__(self):
        return self.__class__(-self.decimal_degree)
    #
    def __pos__(self):
        return self.__class__(self.decimal_degree)
    #
    def __abs__(self):
        self.__pos__()
        return self
    #
    def __add__(self, other, return_class=True):
        """
        Adds something to this dimension
            
            if other is a scalar,
                @return: Another Dimension instance
            if other is another class instance,
                @return: A scalar value 
        """
        if isinstance(other, self.__class__):
            #We're subtracting/adding another class, result should be a scalar
            other = other.decimal_degree
            return_class = False
        elif isinstance(other, BaseDimension):
            #This will start to make sense for derivative classes. If something is an Instance of BaseDimension but not of the same class,
            #Then you've passed in a another type of dimension and should convert it first!
            raise TypeError(u"{my_class} object calculation expected a {my_class}, but you gave a {other_class}. Please convert it to the same type of dimension first.".format(
                            my_class=self.__class__.__name__, other_class=other.__class__.__name__)
                            )
        else:
            other = d(other) #Otherwise, try to cast it to a Decimal
        #
        if return_class: #Assumes we want to return an instance of the class
            return self.__class__(self.decimal_degree + other)
        else:
            return d(self.decimal_degree + other)
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
        Subtracts something from this dimension
            
            if other is a scalar,
                @return: Another Dimension instance
            if other is another class instance,
                @return: A scalar value 
        """
        return self.__add__(-other)
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
        return self.__class__(math.floor(self.decimal_degree))
    #
    def __round__(self):
        return self.__class__(round(self.decimal_degree))
    #   
    def __ceil__(self):
        return self.__class__(math.ceil(self.decimal_degree))
    #   
    def __int__(self):
        return self.degree
    #
    def __float__(self):
        return float(self.decimal_degree)
    #
    def __repr__(self):
        return self.__str__()



class SmartLat(BaseDimension):
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


class SmartLon(BaseDimension):
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


class HourAngle(BaseDimension):
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


class RightAscension(BaseDimension):
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
    

class Altitude(BaseDimension):
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
    

class Azimuth(BaseDimension):
    """
    Represents the Azimuth angle of an AltAz scope. Range 0 to 360,  
    
        0,360  North
        180    South
    """
    base_mode = "dms"
    mutable_mode = False
    str_prefix = "Az"
    range_high = 360
    range_low = 0
    range_rotates = True #True = The range flicks from high back to low if you increase. False = Range oscillates


