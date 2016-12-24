#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot

    Tests
"""
from __future__ import unicode_literals
#
from dateutil import parser
from decimal import Decimal as d
import ephem
import unittest
#
from dimensions import RightAscension, utc_now
from coordinatepairs import ApparentRADec, LatLon
from libraries import sidereal


class TestSingleCoordinates(unittest.TestCase):
    #
    def test_ra_to_ha(self):
        """
        Casts RA to HA and back again, checks they are equal!
        Time is fixed.
        """
        my_time = utc_now()
        my_longitude = "-0.728876"
        target_ra = RightAscension("05","55","10.3053")
        print(u"Right Ascension: %s" % target_ra)
        equiv_ha = target_ra.to_hour_angle(my_longitude, my_time)
        print(u"Equivalent Hour Angle @ %s: %s" % (my_longitude, equiv_ha))
        cast_back_ra = equiv_ha.to_right_ascension(my_longitude, my_time)
        print(u"Casted back Right Ascension: %s" % cast_back_ra)
        #
        self.assertEqual(str(target_ra), str(cast_back_ra))


def compare_pyephem_with_my_code():
    """
    Compares my CoordinatePairs to the Pyephem library
    
    17 Toll Gardens, Bracknell location = 51°24'38.5"N 0°43'44.0"W (51.410684, -0.728876)
    
    Betelgeuse @ 2016-12-23 17:05:
        RA = Decimal('1.5538093862003428125007076232577674090862274169921875')  <radians>
        Dec = Decimal('0.1292627123113075648941361350807710550725460052490234375') <radians>
    """
    latitude = "51:24.385"
    longitude = "-0:43.44"
    my_latlon = LatLon(d("51.410684"), d("-0.728876"))
    my_time = parser.parse("23 Dec 2016 17:15:05")
    #
    # PyEphem's calculations
    observer = ephem.Observer()
    observer.lat = latitude
    observer.lon = longitude
    observer.date = my_time.strftime("%Y/%m/%d %H:%M:%S")
    betelgeuse = ephem.star("Betelgeuse")
    betelgeuse.compute(observer)
    #
    #My calculations - we'll start with Betelgeuse's apparent RADec
    my_radec = ApparentRADec(d(float(betelgeuse.ra)), d(float(betelgeuse.dec)), mode="rad")
    my_radec.apparent = True
    
    #Lets look at output:
    #RA
    print("PyEphem - Betelgeuse RADec = %.6f, %.6f  [%s, %s]" % (betelgeuse.ra, betelgeuse.dec.norm, betelgeuse.ra, betelgeuse.dec))
    print("My - RADec = %.6f, %.6f  %s" % (my_radec.ra.norm, my_radec.dec.norm, my_radec))
    
    #Az alt
    print("PyEphem - Betelgeuse AzAlt = %.6f, %.6f  [%s, %s]" % (betelgeuse.az, betelgeuse.alt, betelgeuse.az, betelgeuse.alt))
    print("My - AzAlt = %s" % my_radec.to_az_alt(latitude=my_latlon.latitude, longitude=my_latlon.longitude, timestamp=my_time))
    
    print("Sidereal RA > %s" % sidereal.raToHourAngle(betelgeuse.ra, my_time, eLong=observer.lon))
    
    print(dir(betelgeuse))
    print(dir(observer))
    


if __name__ == "__main__":
    compare_pyephem_with_my_code()
    #unittest.main()