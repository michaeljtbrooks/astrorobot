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
DEFAULT_ELEVATION = 15 #Metres above sea level

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
EPHEM_OBSERVER_DATE_FORMAT = "%Y/%m/%d %H:%M:%S"
EPHEM_DB_DATE_FORMAT = "%m/%d/%Y %H:%M:%S" #American dates are dumb. Day < month < year... got it? ffs.


#API keys:
KEY_OPEN_WEATHER_MAP = "7b2aaa099558f67a347a178337696d41"
KEY_DARK_SKY = "1611a1fbe455fbed8ba2319280030e3a"