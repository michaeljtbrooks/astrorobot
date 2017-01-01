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
            RADec = a place on the celestial sphere has two modes:
                Actual = the true location if viewed outside Earth's gravity field and outside Earth's atmosphere
                Apparent = where it appears to be when viewed from your location on the Earth
            HaDec = a place on the sky (e.g. a telescope is pointing to) in HourAngle, Declination. This is always relative to the APPARENT RADec.
            AzAlt = a place on the sky (e.g. a telescope is pointing to) in Azimuth, Altitude. This is always relative to the APPARENT RADec. 
                
            

            
            
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
from copy import deepcopy
from datetime import datetime, timedelta
from dateutil import parser
import ephem 
import json
from math import radians
import pytz
import sys
_python3 = sys.version_info > (3,)
if _python3:
    from urllib.parse import urlencode
    from urllib.request import urlopen
else:
    from urllib import urlencode
    from urllib2 import urlopen
#
import settings
from coordinatepairs import RADec, HADec, AzAlt, LatLon, BasePair, ApparentRADec
from dimensions import BaseDimension, SmartLat, ApparentRightAscension, ApparentDeclination
from astrorobot_exceptions import InitError
from utils import d, utc_now
from libraries import sidereal



class TimetrackerMixin(object):
    """
    Carries a time offset from Now. We track time by simply saving an
    offset from true now().
    
    We can freeze this time if we wish by saving the frozen timestamp 
    to a property
    """
    time_travel_offset = None #Gets set to timedelta zero on init
    frozen_time = None #<Datetime> to freeze to
    
    def __init__(self, time_expression=None):
        """
        Creates our time instance:
        @param time_expression: <unicode> or <datetime> The time you wish to start with. If None, assumes now 
        """
        self.timetravel(time_expression)
    
    def __unicode__(self):
        """
        Prints the time
        """
        frozen_flag = ""
        if self.frozen_time:
            frozen_flag = " [⏸ ] FROZEN"
        return "%s%s" % (self.now.strftime("%Y-%m-%d %H:%M:%S"), frozen_flag)
    #
    def __repr__(self, *args, **kwargs):
        return str(self.__unicode__())
    
    def timetravel(self, time=None):
        """
        Time travels you to the specified time!
        
        @param time: The time you wish to emulate.
        
        @return: <datetime> A datetime object for your current time
        """
        now = utc_now()
        if time is None:
            time = now
        elif not isinstance(time, datetime): #Convert whatever has been supplied into a datetime
            time = parser.parse(unicode(time)) #Assume string
        #Make aware and localise:
        if time.tzinfo is None or time.tzinfo.utcoffset(d) is None:
            time = pytz.utc.localize(time)
        #Calculate offset
        self.time_travel_offset = time - now
        return time #Returns the correctly inflated datetime object representing now
    
    def freeze(self, timestamp=None):
        """
        Freezes this time.
        
        @keyword timestamp: If supplied, will freeze the time to this timestamp instead, and continue from it when unfrozen
        
        @return: <datetime> The time we've paused at
        """
        if timestamp: #Update our flag to point to that time
            now = self.timetravel(timestamp)
        else:
            now = self.now
        self.frozen_time = now #This freezes it
        return self.frozen_time
    
    def unfreeze(self):
        """
        Unfreezes the time, carrying on from where the frozen time was
        
        @return: <datetime> The time we've just unpaused from
        """ 
        if self.frozen_time is None: #Do nothing if not frozen
            print("Time on %s() not frozen!!" % (self.__class__.__name__))
            return self.now
        now = self.timetravel(self.frozen_time) #Start the clock ticking again
        self.frozen_time = None #Remove our frozen flag
        return now 
    
    @property
    def now(self):
        """
        Returns the "current" time, tracking seconds elapsed since original time issued
        
        @return: <Datetime> A datetime object for your current time
        """
        if self.frozen_time: #Means we are currently frozen at this time 
            return self.frozen_time
        return utc_now() + self.time_travel_offset
    
    def get_date(self): #Convenience method
        return self.now.date
    
    def get_time(self): #Convenience method
        return self.now.time
    
    @property
    def now_ephem_str(self):
        """
        Gives the time as a string suitable for an ephem Observer to consume
        @return: <str> Representation of the datetime
        """
        out_date_str = self.now.strftime(settings.EPHEM_OBSERVER_DATE_FORMAT)
        return out_date_str
    
    @property
    def now_ephem_db_str(self):
        """
        Gives the time as a string suitable for an ephem database to consume
        @return: <str> Representation of the datetime
        """
        return self.now.strftime(settings.EPHEM_DB_DATE_FORMAT)


class TargetSearchMixin(object):
    """
    Allows us to search for items in the Ephem library and return them
    """
    def ephem_find(self, term, observer=None):
        """
        Searches for an exact ephem item and returns it
        
        @param term: The search term to use
        @keyword observer: The observer to use, falls back to self.ephem_observer 
        
        @return: ephem object, or None if not found
        """
        ephem_target = None
        cap_item = unicode(term).capitalize() #Standardise
        ephem_observer = observer.ephem_observer or self.ephem_observer
        
        #First see if the thing is a solar system object, these are direct properties:
        if cap_item in settings.SOLAR_SYSTEM_OBJECTS:
            ephem_target = getattr(ephem, cap_item)(ephem_observer) # Solar system objects are full classes on ephem
        
        #See if it is a star:
        if not ephem_target:
            try:
                ephem_target = ephem.star(cap_item, ephem_observer)
            except KeyError:
                ephem_target = None
        print("Target: %s" % ephem_target)
        return ephem_target
    #
    def _get_ephem_target(self, item, observer=None, *args, **kwargs):
        """
        Takes a range of inputs and uses it to deduce an ephem_target 
        
        @param target: <Target>, <ephem.Body> or <unicode> - the target we are specifying
        @keyword observer: The observer to use (will ultimately fall back to self.ephem_observer)
        
        @return: ephem.Body - some celestial object
        """
        if isinstance(item, Target): #If a Target is supplied, call its ephem target object
            ephem_target = item.ephem_target
        elif isinstance(item, (ephem.Body, ephem.FixedBody, ephem.EarthSatellite, ephem.EllipticalBody, ephem.Planet, ephem.PlanetMoon, ephem.ParabolicBody)):
            ephem_target = item
        else: #A search term has been supplied, so search the ephem library for it 
            ephem_target = self.ephem_find(item, observer)        
        if ephem_target:
            return ephem_target
        return None


class Observer(TimetrackerMixin, TargetSearchMixin):
    """
    Observer = a person at a position on the earth (Lat/Lon), at a certain time origin (which follows realtime), in specified atmospheric conditions
    
    We use this as a proxy to ephem, allowing us to inject clever adjustments like setting the observer date to TimetrackerMixin.now prior to calling it! 
    """
    ephem_observer = None #Where we'll cache our actual ephem observer
    #
    def __init__(self, latitude=None, longitude=None, location=None, timestamp=None, temperature=settings.DEFAULT_TEMPERATURE, pressure=settings.DEFAULT_PRESSURE, elevation=settings.DEFAULT_ELEVATION, *args, **kwargs):
        """
        Sets up this situation
        
            @keyword latitude: Observer's latitude as a SmartLat object, or as a decimal
            @keyword longitude: Observer's longitude as a SmartLat object, or as a decimal
            -OR-
            @keyword location: Observers LatLon as a LatLon instance, or as a city name
            
            @keyword timestamp: Datetime object reflecting when observing started. If omitted, UTC now will be used
            @keyword temperature: What the current temperature in Celcius of your local area is (for calculating refraction etc)
            @keyword pressure: What the current pressure in mBar of your local area is (for calculating refraction etc) 
        """
        #Sanitise location
        if isinstance(location, (str, unicode)): #User wishes to search for a city:
            try:
                city = ephem.city(location)
            except KeyError: #Location not found in ephem... search online??
                try:
                    city = self.get_ephem_city(location)
                    print("Found online: %s = %s,%s" % (location, city.lat, city.lon))
                except ValueError:
                    raise Exception("Sorry, we cannot find the location '%s' either in our database or online." % location) 
            #Assuming we've found ourselves a little location
            self.location = LatLon(float(city.lat), float(city.lon), mode="rad") #Casts into decimal location
            if city.elevation:
                elevation = d(city.elevation)
        
        if latitude is not None and longitude is not None:
            self.location = LatLon(latitude, longitude) #This will also cast any lat lon strings into proper numbers
        elif location is not None and isinstance(location, LatLon):
            self.location = location
                
        if not location:
            raise InitError(u"You must supply a latitude, longitude or a location (as a LatLon), or valid city name in order to set up an Observer")
        #
        #Sanitise time
        if timestamp is None: #No time supplied, so take now
            self.time_travel_offset = timedelta(seconds=0)
        else:
            #Coerce supplied timestamp into a tz aware datetime object
            #Timestamp has been applied, so convert to offset from now
            self.timetravel(timestamp)
        
        #Resolve the temperature, pressure and elevation:
        self.temperature = temperature or settings.DEFAULT_TEMPERATURE
        self.pressure = pressure or settings.DEFAULT_PRESSURE
        self.elevation = elevation or settings.DEFAULT_ELEVATION
        
        #Create us as an ephem observer:
        self.updated_ephem_observer() #Creates if not already exists
    #
    def __getattr__(self, name, *args):
        """
        This is a proxy method to ensure all attribute calls go through to the 
        underlying ephem_observer object, but with a quick update to its params just before!!
        
        e.g.
            next_rising("sun") = next sun rise
            next_transit("sun") = next sun transit (peak point in sky)
            next_setting("sun") = next sun set
            next_antitransit("sun") = next sun antitransit (lowest point in sky)
            previous_rising("sun") = previous sun rise
            previous_transit("sun") = previous sun transit (peak point in sky)
            previous_setting("sun") = previous sun set
            previous_antitransit("sun") = previous sun antitransit (lowest point in sky)
        
        
        @return: VARIOUS - Depending on the underlying function!
        """
        def wrapper(*args, **kwargs):
            #Quickly update everything (i.e. cast our self properties to ephem.Observer() format):
            ephem_observer = self.updated_ephem_observer()
            
            if unicode(name).startswith("previous_") or unicode(name).startswith("next_"):
                #This is a call to find a target event. Let's go through that first:
                target = self._get_ephem_target(*args, **kwargs) #We search ephem for the ephem_target if applicable
                output = getattr(ephem_observer, name)(target) #We call the function with the target as its first argument
            else:
                #This is a direct proxy through to the ephem.observer, pass in the args and kwargs
                try:
                    output = getattr(ephem_observer, name)(*args, **kwargs) #Assume a function
                except TypeError: #Assume a property
                    output = getattr(ephem_observer, name)
            return output
        return wrapper #This wrapper partial allows us to pass *args and **kwargs into a __getattr__ situation
    #
    def updated_ephem_observer(self):
        """
        Gets our latest ephem observer after it has been updated with the latest model properties
        
        @return: ephem.Observer()
        """
        if self.ephem_observer is None: #Generate a new observer
            self.ephem_observer = ephem.Observer() #Always init with no args
        self.ephem_observer.lat = unicode(self.location.latitude.dd)  #Accepts as a string latitude only
        self.ephem_observer.lon = unicode(self.location.longitude.dd) #Accepts as a string longitude only
        self.ephem_observer.temp = self.temperature
        self.ephem_observer.pressure = self.pressure
        self.ephem_observer.elevation = self.elevation
        self.ephem_observer.date = self.now_ephem_str #Updates the time to the latest tracked value ###THE MOST IMPORTANT BIT!###
        self.ephem_observer.epoch = self.now_ephem_str #Ensures we chase the sky as it looks on today's star atlas, not year 2000
        return self.ephem_observer
    #
    def get_ephem_city(self, address):
        """
        Given a string `address`, do a Google lookup and return an Ephem Observer.
        
        Patched version of ephem.cities.lookup()
        
        @param address: <str> A search string for a place
        
        @return Ephem-observer object with lat, lon, elevation all set
        """
        #First get the location's lat and long
        parameters = urlencode({'address': address, 'sensor': 'false'})
        url = 'http://maps.googleapis.com/maps/api/geocode/json?' + parameters
        data = json.loads(urlopen(url).read().decode('utf-8'))
        results = data['results']
        if not results:
            raise ValueError('Google cannot find a place named %r' % address)
        address_components = results[0]['address_components']
        location = results[0]['geometry']['location']
        
        #Now fetch the elevation:
        parameters = urlencode({'locations': "%s,%s" % (location["lat"], location["lng"])})
        url = 'http://maps.googleapis.com/maps/api/elevation/json?' + parameters
        data = json.loads(urlopen(url).read().decode('utf-8'))
        print(data)
        elev_results = data['results']
        elevation = settings.DEFAULT_ELEVATION
        if elev_results:
            elevation = elev_results[0]['elevation']
        
        #Create our ephem_observer
        o = ephem.Observer()
        o.name = ', '.join(c['long_name'] for c in address_components)
        o.lat = radians(location['lat'])
        o.lon = radians(location['lng'])
        o.elevation = elevation
        return o
    #
    def get_weather(self, timestamp=None):
        """
        This fetches the weather for your location, at the given time (defaults to now).
        Future dates will give a forecast. Past dates will give nearest actual observations. 
        
        @keyword timestamp: <datetime> The time you wish to explore. Defaults to self.now
        
        @return: {
                    temperature : <temperature Celcius>,
                    pressure : <pressure millibar>,
                    humidity : <humidity %>,
                    cloud_cover : <cover %>,
                    transparency : <?units>,
                    seeing : <int> 1-5 Antoniadi scale,
                }
        
        @TODO:
        """
        pass
    #
    def _get_sidereal_time_obj(self, timestamp=None, longitude=None):
        """
        Resolves the sidereal time of the specified location at the specified timestamp
        
        @keyword timestamp: <datetime> The UTC time you wish to query for, defaults to self.now (Observer's time)
        @keyword longitude: <LatLon> or <SmartLon> The longitude you wish to query for (defaults to self.location.longitude) 
        
        @return: <SiderealTime> object representing your local hour angle
            
            (NB: SiderealTime instance carries two useful properties:
                .hours: the hours as a float
                .radians: the hours in radians (needs normalising to 2*Pi
        
        """
        if timestamp is None:
            timestamp = self.now
        if longitude is None:
            longitude_radians = self.location.longitude.radians()
        
        greenwich_sidereal_time = sidereal.SiderealTime.fromDatetime(timestamp) #The Sidereal Time at longitude 0.0
        location_sidereal_time = greenwich_sidereal_time.lst(longitude_radians) #Your location's Sidereal Time (at your longitude)
        
        return location_sidereal_time
    #
    def get_sidereal_time(self, timestamp=None, longitude=None):
        """
        Resolves the sidereal time of the specified location at the specified timestamp
        
        @keyword timestamp: <datetime> The UTC time you wish to query for, defaults to self.now (Observer's time)
        @keyword longitude: <LatLon> or <SmartLon> The longitude you wish to query for (defaults to self.location.longitude) 
        
        @return: <Decimal> sidereal time in hours
        """
        lst_time_out = self._get_sidereal_time_obj(timestamp=timestamp, longitude=longitude).hours #Returns a float hours
        return d(lst_time_out) #Return as Decimal
    #
    @property
    def st(self):
        return self.get_sidereal_time()
    @property
    def sidereal(self):
        return self.get_sidereal_time()
    @property
    def sidereal_time(self):
        return self.get_sidereal_time()



class Target(TargetSearchMixin):
    """
    Represents an object on the celestial sphere, may be star / planet / satellite / piece of space junk
    
    Wrapper around on PyEphem's skyobject types.
    
    Will contain measurements for:
        a_ra = astrometric RightAscension
        a_dec = astrometric Declination
        a_ra_dec = astrometric RADec() as a pair
        ra = APPARENT RightAscension
        dec = APPARENT Declination
        ra_dec = APPARENT RADec() as a pair
        
        ha = APPARENT HourAngle
        dec = APPARENT Declination
        ha_dec= APPARENT HADec as a pair
        
        az = APPARENT azimuth
        alt = APPARENT altitude
        az_alt = APPARENT AzAlt as a pair
        
    """
    ephem_target = None #Will store PyEphem's object here
    observer = None #Stores our observer
    ephem_observer = None #Stores our Ephem Observer object, allowing us to use common mixin
    
    def __init__(self, observer=None, name=None, ra=None, dec=None, az=None, alt=None, ha=None, target_type="fixed"):
        """
        TARGET: Creates the relevant target either by lookup name or coordinates
        
        @TODO: build in creating new object from coords
        
        REQUIRED:
        @param observer: <Observer> An observer (person on a place on the earth at a certain time)
        
        TO SPECIFY THE TARGET's POSITION
        @keyword name:  <unicode> The name of the item, which will be used to search the ephem database
                        <Target> An already found Target object which can be used with a new observer
        -OR-
        @keyword ra: <RightAscension> or (tuple: hr,min,sec) or <Decimal(RA)> The target's position in uncorrected astrometric RightAscension        
        @keyword dec: <Declination> or (tuple: deg,min,sec) or <Decimal(dec)> The target's position in uncorrected astrometric Declination
        -OR-
        @keyword ha: <HourAngle> or (tuple: hr,min,sec) or <Decimal(HA)> The target's position in uncorrected astrometric HourAngle        
        @keyword dec: <Declination> or (tuple: deg,min,sec) or <Decimal(dec)> The target's position in uncorrected astrometric Declination
        -OR-
        @keyword az: <Azimuth> or (tuple: deg,min,sec) or <Decimal(Az)> The target's position in azimuth        
        @keyword alt: <Altitude> or (tuple: deg,min,sec) or <Decimal(Alt)> The target's position in altitude
        """
        if not isinstance(observer, Observer):
            raise InitError("Target must be passed an Observer to know where it is!")
        self.observer = observer
        self.ephem_observer = observer.ephem_observer
        
        print("Target init name: %s " % name)
        if name:
            ephem_target = self.ephem_find(name, observer)
            self.ephem_target = ephem_target
        
    def __getattr__(self, name, *args):
        """
        PROXY to ephem.Body
        
        This is a proxy method to ensure all attribute calls go through to the 
        underlying ephem_Body object, but with a quick update to the observer's time just before
        
        @return: VARIOUS - Depending on the underlying function!
        """
        self.ephem_observer = self.observer.updated_ephem_observer() #Update the latest time
        self.ephem_target.compute(self.ephem_observer) #Updates this to the latest observer time
        output_attr = getattr(self.ephem_target, name)
        if callable(output_attr): #Proxy to a wrapper function if this is a callable
            def wrapper(*args, **kwargs):
                output = getattr(self.ephem_target, name)(*args, **kwargs)
                return output
            return wrapper #This wrapper partial allows us to pass *args and **kwargs into a __getattr__ situation
        else: #If not a function, just return the attr
            return output_attr
    
    def position_at(self, timestamp):
        """
        Returns this target's expected position at the specified time
        
        @param timestamp: <datetime> UTC time when you'd like to determine the object's position for 
        
        @return: <Target> updated to point to specified time
        """
        new_target = deepcopy(self) #Make a copy so we don't pollute our original
        new_target.observer.timetravel(timestamp) #This will be a new observer object too, yes things can get complicated if an observer moves!!
        return new_target #Any properties called on that daughter object will 
    
    #--- Dimensions and CoordinatePairs ---
    @property
    def ra_dec(self):
        """
        Returns the RightAscension & Declination as a Coordinate pair
        
        @return <AzAlt> of target
        """
        return ApparentRADec(self.ra, self.dec, mode="rad")
    
    @property
    def hour_angle(self):
        """
        Converts this target's Apparent RA to Hour Angle for the position of the observer
        
        @return: <HourAngle> of target 
        """
        target_apparent_ra = ApparentRightAscension(self.ra, mode="rad")
        return target_apparent_ra.to_hour_angle(self.observer.location.longitude, self.observer.now)
    @property
    def ha(self):
        return self.hour_angle

    @property
    def ha_dec(self):
        """
        Returns the HourAngle, Declination pair for this target. Suitable for pointing a scope to!
        
        @return: <HADec> of target
        """
        target_ha = self.hour_angle
        target_apparent_dec = ApparentDeclination(self.dec, mode="rad")
        return HADec(target_ha, target_apparent_dec)
    
    @property
    def az_alt(self):
        """
        Returns the Azimuth & Altitude as a Coordinate pair
        
        @return <AzAlt> of target
        """
        return AzAlt(self.ephem_target.az, self.ephem_target.alt, mode="rad")



class BaseSituation(object):
    """
    ABSTRACT CLASS
    
    Holds the relevant components for a situation = an observer on the Earth using a device to look at the celestial sphere
    This class brings together the following ingredients:
        observer = Observer (a person starting to look at the sky from a certain time at a certain geographical location)
        pointing_target = Target (currently being pointed to)
        sought_target = Target (being sought)
        telescope = TelescopeWithMount (represents the telescope and the mount)
            focus
            motors
            filter
            lenses
        camera = Camera (the optical acquisition device being used)
    
     
    
    equatorial = HADec instance
    celestial = RADec instance
    apparent_celestial = RADec instance corrected for refraction etc
    azalt = AzAlt instance
    location = Position on the Earth where the observer is
    time_travel_offset = Timedelta from current real UTC that this situation reflects
    observer = PyEphem observer instance which is used internally to compute apparent HAdec etc
    """
    observer = None #PyEphem observer
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
        observer.date = timestamp.strftime(settings.EPHEM_OBSERVER_DATE_FORMAT)
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
        self.observer.date = time_now.strftime(settings.EPHEM_OBSERVER_DATE_FORMAT) #Make sure the observer knows the time has changed
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
        observer.date = time_now.strftime(settings.EPHEM_OBSERVER_DATE_FORMAT)
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
            self.observer.date = dt_now.strftime(settings.EPHEM_OBSERVER_DATE_FORMAT)
        
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
        cap_item = unicode(item).capitalize() #Standardise
        
        observer = deepcopy(self.observer) #Copy our observer to avoid polluting self.observer
        
        if timestamp is None:
            time_now = self.now
        else:
            time_now = timestamp #We do not update self.time_travel_offset here
        
        #First see if the thing is a solar system object, these are direct properties:
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
        astronomic_ra_dec = RADec(target.ra.split(":"), target.dec.split(":"))
        apparent_ra_dec = astronomic_ra_dec.apparent()
        return apparent_ra_dec
  
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
        pass


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
