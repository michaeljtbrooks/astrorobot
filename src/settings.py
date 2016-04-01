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