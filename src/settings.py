#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    Settings
"""
from decimal import Decimal as d
#from pint import UnitRegistry as U #For expressing physical quantities

### User-specific constants ###
LOCATION_LAT_LON = (51.410684, -0.728876)

### Scientific Constants ###
EARTH_ROTATION_SIDEREAL = 1

### Program defaults ###
DEFAULT_TEMPERATURE = 8.0 #Celsius (as a float)
DEFAULT_PRESSURE = 1010 #mBar (as int)
DEFAULT_ELEVATION = 0 #Metres above sea level

## Built-in celestial objects ###
SOLAR_SYSTEM_OBJECTS = (
    'Ariel',
    'B1900',
    'B1950',
    'Callisto',
    'Deimos',
    'Dione',
    'EarthSatellite',
    'Enceladus',
    'Europa',
    'Ganymede',
    'Hyperion',
    'Iapetus',
    'Io',
    'Jupiter',
    'Mars',
    'Mercury',
    'Mimas',
    'Miranda', 
    'Moon',
    'Neptune', 
    'Oberon',
    'Observer',
    'Phobos',
    'Pluto',
    'Rhea',
    'Saturn',
    'Sun',
    'Tethys',
    'Titan',
    'Titania',
    'Umbriel',
    'Uranus',
    'Venus',
)

## How our date is formatted
EPHEM_OBSERVER_DATE_FORMAT = "%Y/%m/%d %H:%i:%s"
EPHEM_DB_DATE_FORMAT = "%m/%d/%Y %H:%i:%s" #American dates are dumb. Day < month < year... got it? ffs.