#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    
    Reorganisation:
        - Change BaseSituation to "Observer", be a wrapper around PyEphem's Observer. Time can be a property of the Observer
        - Write "Target" class = PyEphem's Target, but with coordinates cast into our native types
        - Create classes for the interchangeable components of the telescope setup
            - Camera
            - Motors
            - Mount & Scope
        - Create "Equipment" class which brings together the camera, motors, mount & scope (?reads from ext config file)
        - Create  proper "Situation" class which brings together:
            - Equipment: (Mounts & Scope + motors + camera, inc whether position has been declared)
            - Observer
            - Target (what you're currently pointing at or tracking, if known)
    
    
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
                uncorrected_celestial = (RA, Declination)                   #Where on the skymap your mount is currently pointing (varies unless tracking)
                apparent_celestial = (RA, Declination)          #Where on the skymap your mount appears to be pointing when refraction taken into account (varies unless tracking)
                azalt = (azimuth, altitude)                     #Where an azalt mount is currently pointing (fixed unless motors move)
                location = (latitude, longitude)                #The observer's position on the Earth
                time_travel_offset = (timestamp)                #How far forwards or back from UTC now is the observer
            
                prime_subjective = Which axis we are using to calculate all others from
                humidity = The current local humidity (affects the refraction cooefficient which therefore affects apparent_celestial position)
                

"""
from __future__ import unicode_literals
#
from datetime import datetime, timedelta
from dateutil import parser
import ephem
import pytz
#
from coordinatepairs import RADec, HADec, AzAlt, LatLon, BasePair
from dimensions import BaseDimension
import settings
from utils import d, utc_now




class InitError(SyntaxError):
    """
    A class to represent incorrect setting up of your situation 
    """
    pass


class Target(object):
    """
    Represents an object on the celestial sphere, may be star / planet / satellite / piece of space junk
    
    Wrapper around on PyEphem's skyobject types,
    """
    sky_object = None #Will store PyEphem's object here
    pass





class BaseSituation(object):
    """
    TODO: Change this to "Observer", be a wrapper around PyEphem's Observer
    
    ABSTRACT CLASS
    
    Holds the relevant components for a situation = an observer on the Earth using a device to look at the celestial sphere
    Wraps around the PyEphem Observer property, but will automatically update the observer's date to "now" such that time is tracked 
    
    equatorial = HADec instance
    celestial = RADec instance
    apparent_celestial = RADec instance corrected for refraction etc
    azalt = AzAlt instance
    location = Position on the Earth where the observer is
    time_travel_offset = Timedelta from current real UTC that this situation reflects
    observer = PyEphem observer instance which is used internally to compute apparent HAdec etc
    """
    ephem_observer = None
    prime_subjective = "ha_dec" #Which element to calculate the others from / to affect if move() called
    initialised = False #Marks when this situation has been started
    #
    def __init__(self, latitude=None, longitude=None, location=None, timestamp=None, pointing=None, temperature=settings.DEFAULT_TEMPERATURE, pressure=settings.DEFAULT_PRESSURE, *args, **kwargs):
        """
        Sets up this situation
        
            @keyword latitude: Observer's latitude as a SmartLat object, or as a decimal
            @keyword longitude: Observer's longitude as a SmartLat object, or as a decimal
            -OR-
            @keyword location: Observers LatLon as a LatLon instance
            
            @keyword timestamp: Datetime object reflecting when observing started. If omitted, UTC now will be used
            
            @keyword pointing: What your scope is pointing to as a RADec or HAdec or AzAlt instance. Default = 0,0
            
            @keyword temperature: What the current temperature in Celcius of your local area is (for calculating refraction etc)
            
            @keyword pressure: What the current pressure in mBar of your local area is (for calculating refraction etc) 
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
        #From now on, get time as self.now
        
        #Resolve the temperature, pressure and elevation:
        self.temperature = temperature or settings.DEFAULT_TEMPERATURE
        self.pressure = pressure or settings.DEFAULT_PRESSURE
        self.elevation = settings.DEFAULT_ELEVATION
        
        #Resolve pointing into the three vars:
        if not isinstance(pointing, BasePair):
            print(u"WARNING: Could not make sense of the direction you supplied for where you are pointing the telescope. Assuming at Zenith.")
            pointing = AzAlt(0,90)
        self.az_alt = pointing.to_az_alt(latitude, longitude, timestamp)
        self.ha_dec = pointing.to_ha_dec(latitude, longitude, timestamp) #THIS quantity should adjust to apparent!!!
        self.ra_dec = pointing.to_ra_dec(latitude, longitude, timestamp)
        
        #Build our observer:
        observer = self.build_pyephem_observer(location, temperature, pressure, timestamp)
        
        self.ra_dec_apparent = self.ra_dec.apparent(observer=observer, latitude=latitude, longitude=longitude, timestamp=timestamp, temperature=temperature, pressure=pressure) 
        #
        #Now mark initialisation complete:
        self.initialised = True
    #
    def build_pyephem_observer(self, location, temperature=None, pressure=None, timestamp=None):
        """
        Creates a pyephem Observer instance, which we can use to compute apparent RAdec thus HAdec
        
        Uses functions from PyEphem [http://rhodesmill.org/pyephem/index.html]
        
        @param location: A LatLon instance of the observer
        @keyword temperature: Observing temperature. If omitted will use default
        @keyword pressure: Observing pressure. If omitted, will use default
        @keyword timestamp: Observing time. If omitted, will default to "now" (i.e. observer's time, ticking away)
        
        @returns PyEphem observer instance
        """
        if timestamp is None:
            timestamp = self.now
        observer = ephem.Observer()
        observer.lat = unicode(location.latitude.dd) #Accepts as a string only
        observer.lon = unicode(location.longitude.dd)#Accepts as a string only
        observer.date = timestamp.strftime(settings.EPHEM_DATE_FORMAT)
        observer.temp = temperature
        observer.pressure = pressure
        self.observer = observer
        return observer
    
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
        self.observer.date = time_now.strftime(settings.EPHEM_DATE_FORMAT) #Make sure the observer knows the time has changed
        self.az_alt = master_pointing.to_az_alt(self.latitude, self.longitude, time_now)
        self.ha_dec = master_pointing.to_ha_dec(self.latitude, self.longitude, time_now)
        self.ra_dec = master_pointing.to_ra_dec(self.latitude, self.longitude, time_now)
        self.ra_dec_apparent = self.ra_dec.apparent(observer=self.observer, latitude=self.latitude, longitude=self.longitude, timestamp=time_now, temperature=self.temperature, pressure=self.pressure)
        return self
    #
    def update_observer_time(self, observer=None, timestamp=None):
        """
        Update the date and time of the observer. Defaults to self.now
        
        @keyword observer: The observer, If omitted, will use self.observer
        @keyword timestamp: The time to move the observer to. If omitted, will use self.now
        
        @return: <Observer> 
        """
        observer = observer or self.observer
        if timestamp:
            time_now = self.timetravel(timestamp)
        else:
            time_now = self.now
        observer.date = time_now.strftime(settings.EPHEM_DATE_FORMAT)
        return observer
    #
    def change_observer_position(self, location=None, latitude=None, longitude=None, timestamp=None, temperature=None, pressure=None, elevation=None):
        """
        Changes the position of an observer to the provided values:
        
        @keyword location: <LatLon> GPS coordinates
        @keyword latitude: decimal latitude. If None will not be updated
        @keyword longitude: decimal longitude. If None will not be updated
        @keyword timestamp: <Datetime> timetravel to the appropriate point. If None will be ignored
        @keyword temperature: <decimal> temperature in Celsius. If None will be ignored
        @keyword elevation: <decimal> metres ASL. If None, will be ignored
        
        @return: self  
        """
        #Deal with position:
        if latitude is not None and longitude is not None:
            self.location = LatLon(latitude, longitude)
            self.observer.lat = latitude
            self.observer.lon = longitude
        elif location is not None and isinstance(location, LatLon):
            self.location = location
            self.observer.lat = location.latitude
            self.observer.lon = location.longitude
        
        #Timetravel to the time provided:
        if timestamp:
            dt_now = self.timetravel(timestamp)
            self.observer.date = dt_now.strftime(settings.EPHEM_DATE_FORMAT)
        
        #Deal with other params:
        if pressure is not None:
            self.observer.pressure = pressure
        if temperature is not None:
            self.observer.temp = temperature
            
        #Elevation not yet implemented
        
        self.update() #Adjust our non-primary vals to suit the new location
        
        return self
    #
    def set_pointing_to(self, x=None, y=None):
        """
        Moves the pointing position by the given number of decimal degree units, or to the specified target
        You can either pass in a basepair as the first value, or dimensions / scalars as two values
        
        @param x:  <BasePair> | <BaseDimension> | <decimal> scalar value
        @keyword y: <BaseDimension> | <decimal> scalar value
        
        @return: self
        
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
            return self.set_pointing_to(target_pointing)
        else:
            #User has probably passed in two scalar values, they will want to MOVE the scope by the amount
            if x is None:
                x = d(0)
            else: #Force the x into decimal
                x = d(x)
            if y is None:
                y = d(0)
            else:
                y = d(y) #Force into decimal
            current_pointing = self.get_primary()
            new_pointing = current_pointing + (x,y) #All BasePairs should be able to deal with a tuple of scalars, returning an instance of their class again
            self.set_primary(new_pointing) #Set our primary value to this
        self.update() #Updates all secondary (dependent) measures
        return self #For chaining 
    #
    def where_is(self, item, timestamp=None, mode=None):
        """
        Tells you where a certain celestial object APPARENTLY is *at self.now*
        
        @param item: <unicode> A description of the item you are seeking
        @keyword timestamp: <datetime> when we are looking. Defaults to time now 
        @keyword mode: What to cast the output coords into (RADec / HA / AltAz)
        
        @return: <RADec> Apparent coordinates of the item
        """
        cap_item = unicode(item).capitalize()
        
        if timestamp is None:
            timestamp = self.now
        
        #First see if the thing is a solar system object:
        if cap_item in settings.SOLAR_SYSTEM_OBJECTS:
            target = getattr(ephem, cap_item)(self.observer, epoch=time_now) # Solar system objects are full classes on ephem
        
        #See if it is a star:
        if not target:
            try:
                target = ephem.star(cap_item, self.observer)
            except KeyError:
                target = None
        
        if not target: #Means item not found
            return None
        
        #Ask ephem to compute the APPARENT ra_dec
        target.compute(time_now, epoch=time_now)
        uncorrected_ra_dec = RADec(target.ra.split(":"), target.dec.split(":"))
        corrected_ra_dec = uncorrected_ra_dec.apparent()
        
        
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
    def apparent_celestial(self):
        return self.ra_dec_apparent
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
