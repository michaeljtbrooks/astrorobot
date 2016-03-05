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
        
            A CelestialPosition has three components:
                equatorial = (hour_angle,Declination)           #Where to point your equatorial mount (varies with time)
                    - hour_angle is in HOURS
                celestial = (RA, Declination)                   #Where object is on the skymap (fixed)
                    - RA is in HOURS
                azalt = (azimuth, altitude)                     #Where to point your AzAlt mount (varies with time)
                    - azimuth is in DEGREES
            
            An ObservingPosition has five components:
                equatorial = (hour_angle,Declination)           #Where the equatorial mount is currently pointing (fixed unless motors move)
                celestial = (RA, Declination)                   #Where on the skymap your mount is currently pointing (varies unless tracking)
                azalt = (azimuth, altitude)                     #Where an azalt mount is currently pointing (fixed unless motors move)
                location = (latitude, longitude)                #The observer's position on the Earth
                time = (timestamp)                              #The UTC time of observation

"""

from coordinatepairs import RADec, HADec, AzAlt

class BaseSituation(object):
    """
    ABSTRACT CLASS
    
    Holds the relevant components for a situation
    
    equatorial = HADec instance
    celestial = RADec instance
    azalt = AzAlt instance
    """
    variables = None
    fixed = ["ra_dec", "ha_dec", "az_alt"]
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


class CelestialPosition(BaseSituation):
    """
    A position on the Celestial Sphere. Invariant in time, longitude and latitude,
    
        A CelestialPosition has three components:
            equatorial = (hour_angle,Declination)           #Where to point your equatorial mount (varies with time)
                - hour_angle is in HOURS
            celestial = (RA, Declination)                   #Where object is on the skymap (fixed)
                - RA is in HOURS
            azalt = (azimuth, altitude)                     #Where to point your AzAlt mount (varies with time)
                - azimuth is in DEGREES
    """
    def __init__(self, *args, **kwargs):
        """
        Set up by passing in an RADec object, or RA,Dec in degrees
        
        @param: An RADec instance or derivative class instance
        -OR-
        @param x: A RA dimension class or degree value
        @param y: A Dec dimension class or degree value
        """
        if isinstance(args[0], RADec): #Already fully inflated object there
            self.ra_dec = args[0]
        else:
            self.ra_dec = RADec(args[0], args[1])
    #
    def update_by_RA(self, latitude, longitude, timestamp=None):
        """
        Takes RA as gospel, then changes the other two coordinates to match 
        """
        ##HERE##
        pass


class ObservingPosition(BaseSituation):
    pass