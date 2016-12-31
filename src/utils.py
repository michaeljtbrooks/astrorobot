#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    Utils
        Provides useful methods to use across location / coordinate objects
        
    
    Nomenclature            
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
from datetime import datetime
from decimal import Decimal 
import math
#
import pytz
from ephem import Angle


### Constants ###
PI = Decimal("3.141592653589793238462643383279502884197169399") #I hate floats! This is a far more accurate Decimal Pi
TWO_PI = Decimal("2") * PI
TWO_PI_FLOAT = float(TWO_PI)


def d(val, *args, **kwargs):
    """
    Casts the value into a true decimal
    """
    if not isinstance(val, Decimal):
        val_str = unicode(val)
        val_d = Decimal(val_str, *args, **kwargs)
    else:
        val_d = val
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
    if isinstance(radians, (Angle,)): #Ephem doesn't output a numerical value with unicode
        radians = Decimal(float(radians))
    else:
        radians = Decimal(unicode(radians))
    degrees = radians * (360/TWO_PI)
    degrees = degrees % (Decimal("360"))
    return degrees

def coord_rotate_rad(x, y, z):
    """Used to convert between equatorial and horizon py.
        
        Adapted from http://www.nmt.edu/tcc/help/lang/python/examples/sidereal/ims/
        
        For Converting Az,Alt > RA,Dec:
            x = Altitude in radians (-pi to pi)
            y = Latitude in radians (-pi to pi)
            z = Azimuth in radians (0 to 2pi)
            
            > xt = Declination in radians
            > yt = HourAngle in radians
            
        For Converting RA,Dec > Az,Alt:
            x = Declination in radians (-pi to pi)
            y = Latitude in radians (-pi to pi)
            z = HourAngle in radians (0 to 2pi), (which requires longitude if converting from RA)
            
            > xt = Altitude in radians
            > yt = Azimuth in radians
        
      [ x, y, and z are angles in radians ->
          return (xt, yt) where
          xt=arcsin(sin(x)*sin(y)+cos(x)*cos(y)*cos(z)) and
          yt=arccos((sin(x)-sin(y)*sin(xt))/(cos(y)*cos(xt))) ]
    """
    #-- 1 --
    xt  =  math.asin ( math.sin(x) * math.sin(y) +
                  math.cos(x) * math.cos(y) * math.cos(z) )
    #-- 2 --
    yt  =  math.acos ( ( math.sin(x) - math.sin(y) * math.sin(xt) ) /
                  ( math.cos(y) * math.cos(xt) ) )
    #-- 3 --
    if  math.sin(z) > 0.0:
        yt  =  TWO_PI - yt

    #-- 4 --
    return (xt, yt)