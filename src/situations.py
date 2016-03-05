#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    SITUATIONS
        
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

from datetime import datetime, timedelta
from dateutil import parser
import pytz
#
from coordinatepairs import RADec, HADec, AzAlt, LatLon, BasePair
from dimensions import BaseDimension
from utils import d, utc_now




class InitError(SyntaxError):
    """
    A class to represent incorrect setting up of your situation 
    """
    pass


class BaseSituation(object):
    """
    ABSTRACT CLASS
    
    Holds the relevant components for a situation = an observer on the Earth using a device to look at the celestial sphere
    
    equatorial = HADec instance
    celestial = RADec instance
    azalt = AzAlt instance
    location = Position on the Earth where the observer is
    time_travel_offset = Timedelta from current real UTC that this situation reflects
    """
    prime_subjective = "ha_dec" #Which element to calculate the others from / to affect if move() called
    initialised = False #Marks when this situation has been started
    #
    def __init__(self, latitude=None, longitude=None, location=None, timestamp=None, pointing=None, *args, **kwargs):
        """
        Sets up this situation
        
            @keyword latitude: Observer's latitude as a SmartLat object, or as a decimal
            @keyword longitude: Observer's longitude as a SmartLat object, or as a decimal
            -OR-
            @keyword location: Observers LatLon as a LatLon instance
            
            @keyword timestamp: Datetime object reflecting when observing started. If omitted, UTC now will be used
            
            @keyword pointing: What your scope is pointing to as a RADec or HAdec or AzAlt instance. Default = 0,0
        """
        #Sanitise location
        if latitude is not None and longitude is not None:
            self.location = LatLon(latitude, longitude)
        elif location is not None and isinstance(location, LatLon):
            self.location = location
        else:
            raise InitError(u"You must supply a latitude, longitude or a location (as a LatLon) in order to set up an astronomy Situation")
        #
        #Sanitise time
        now = utc_now()
        if timestamp is None: #No time supplied, so take now
            self.time_travel_offset = timedelta(seconds=0)
        else:
            #Coerce supplied timestamp into a tz aware datetime object
            #Timestamp has been applied, so convert to offset from now
            self.timetravel(timestamp)
        #
        #Resolve pointing into the three vars:
        if not isinstance(pointing, BasePair):
            print(u"WARNING: Could not make sense of the direction you supplied for where you are pointing the telescope. Assuming at Zenith.")
            pointing = AzAlt(0,90)
        self.az_alt = pointing.to_az_alt(latitude, longitude, timestamp)
        self.ha_dec = pointing.to_ha_dec(latitude, longitude, timestamp)
        self.ra_dec = pointing.to_ra_dec(latitude, longitude, timestamp)
        #
        #Now mark initialisation complete:
        self.initialised = True
    #
    def get_primary(self):
        """
        Returns the primary measure instance
        """
        return getattr(self, self.prime_subjective) #Our prime subjective measure (e.g. self.az_alt)
    #
    def set_primary(self, value):
        """
        Sets the primary measure instance to something
        """
        setattr(self, self.prime_subjective, value)
        return self
    #
    def update(self):
        """
        Updates all the secondary measures according to the primary measure
        
        @return: self (for chaining)
        """
        master_pointing = self.get_primary()
        time_now = self.now #Means "now" according to this observer's time
        self.az_alt = master_pointing.to_az_alt(self.latitude, self.longitude, time_now)
        self.ha_dec = master_pointing.to_ha_dec(self.latitude, self.longitude, time_now)
        self.ra_dec = master_pointing.to_ra_dec(self.latitude, self.longitude, time_now)
    #
    def move(self, x=None, y=None):
        """
        Moves the pointing position by the given number of decimal degree units, or to the specified target
        """
        if isinstance(x, BasePair):
            #User has supplied a Coordinate pair to move to. Need to cast it into our correct position
            cast_to_prime_str = "to_%s" % self.prime_subjective #This is what the method is to cast to our target dimension
            new_pointing = getattr(x, cast_to_prime_str)(self.latitude, self.longitude, self.now) #Turns the input basepair into the correct format
            self.set_primary(new_pointing) #Set our primary measure to this BasePair
        elif isinstance(x, BaseDimension):
            #User has passed in a pair of dimensions. We need to find out what these are, then combine them to get a point, then call move() again
            try:
                x_type = x.str_prefix
                y_type = y.str_prefix
            except (TypeError,AttributeError):
                raise TypeError(u"Situation.move(X,Y): You passed in a Dimension instance for X, but something else ({y_class}) for Y. Please enter either a CoordinatePair, or a compatible set of Dimensions, or scalar values to move().".format(y_class=y.__class__.__name__))
            if x_type == "HrA" and y_type == "Dec": #It's an HADec pair
                target_pointing = HADec(x,y)
            elif x_type == "RA" and y_type == "Dec": #It's an RADec pair
                target_pointing = RADec(x,y)
            elif x_type == "Az" and y_type == "Alt": #It's an AzAlt pair
                target_pointing = AzAlt(x,y)
            else:
                #Incompatible pair
                raise TypeError(u"Situation.move(X,Y): You passed in incompatible dimensions for X and Y. Please enter either a CoordinatePair, or a compatible set of Dimensions, or scalar values to move().")
            #Now call move again on our target_pointing (this time the method will process it as a pair:
            return self.move(target_pointing)
        else:
            #User has probably passed in two scalar values, they will want to MOVE the scope by the amount
            if x is None:
                x = d(0)
            if y is None:
                y = d(0)
            current_pointing = self.get_primary()
            new_pointing = current_pointing + (x,y) #All BasePairs should be able to deal with a tuple of scalars, returning an instance of their class again
            self.set_primary(new_pointing) #Set our primary value to this
        self.update() #Updates all secondary (dependent) measures
        return self #For chaining 
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
        self.time_travel_offset = time - utc_now()
        #Update other values if initialised
        if self.initialised:
            self.update()
        return time #Returns the correctly inflated datetime object representing now    
    #
    #PROPERTY ALIASES
    @property
    def equatorial(self):
        return self.ha_dec
    #
    @property
    def celestial(self):
        return self.ra_dec
    #
    @property
    def azalt(self):
        return self.az_alt
    #
    @property
    def latitude(self):
        return self.location.latitude
    #
    @property
    def longitude(self):
        return self.location.longitude
    #
    @property
    def now(self):
        """
        Returns the "current" observer's time once epoch offset taken into account
        """
        return utc_now() + self.time_travel_offset


class ObservingPositionEquatorial(BaseSituation):
    """
    An observer using an equatorially mounted scope
    """
    prime_subjective="ha_dec"


class ObservingPositionHorizontal(BaseSituation):
    """
    An observer using an AzAlt type mount
    """
    prime_subjective="az_alt"
