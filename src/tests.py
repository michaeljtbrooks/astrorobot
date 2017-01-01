#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot

    Tests
"""
from __future__ import unicode_literals
#
from dateutil import parser
from datetime import datetime
from decimal import Decimal as d
import ephem
import unittest
#
from dimensions import RightAscension, utc_now
from coordinatepairs import ApparentRADec, LatLon
from libraries import sidereal
from situations import Observer, Target
import time
import settings


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
    my_time = datetime.utcnow()
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


def give_it_a_spin():
    """
    Tests basic functionality of Astrorobot
    
    See: https://github.com/Stellarium/stellarium/blob/585dcacdf71372de845b81420e1b9d1b80369263/src/core/RefractionExtinction.cpp
        for calculations on
            Extinction = amount of light absorbed by atmosphere vs altitude (NB also causes redenning of objects as blue light absorbed more)
            Refraction = bending of light by atmosphere vs altitude (more bending nearer horizon
            
    This is the refraction calculating code used by libastro (which ephem uses):
        https://fossies.org/dox/xephem-3.7.7/refract_8c_source.html
    """
    dt = datetime(2016,12,31,17,25,00)
    ob = Observer(location="Bracknell,UK", timestamp=dt)
    s = Target(ob, "Betelgeuse")
    print("Betelgeuse:")
    print("\tEphem a-RADec: %s, Dec: %s  (expected: 5:56:6 / 7:24:32)" % (s.a_ra, s.a_dec))
    print("\tEphem RADec: %s, Dec: %s" % (s.ra, s.dec))
    print("\tEphem AzAlt: Az: %s, Alt: %s" % (s.az, s.alt))
    print("\tAstrorobot RADec apparent: %s" % unicode(s.ra_dec))
    print("\tAstrorobot HADec apparent: %s  (expected: 18:08:28 / 7:30:10)" % unicode(s.ha_dec))
    print("\tAstrorobot AzAlt apparent: %s  (expected: 86:56:56 / 7:10:32)" % unicode(s.az_alt))
    #And now check directly
    e_o = ephem.Observer()
    e_o.lat = 51.410684
    e_o.lon = -0.728876
    e_o.elev = 50
    e_o.horizon = 0
    e_o.date = dt.strftime(settings.EPHEM_OBSERVER_DATE_FORMAT)
    e_o.pressure = 0 #Set to zero
    e_o.temp = 10.0
    t = ephem.star("Betelgeuse")
    t.compute(e_o)
    print("Ephem pressure 0 AzAlt: %.12f, %.12f  |  RADec: %.12f, %.12f" % (t.az, t.alt, t.ra, t.dec))
    e_o.pressure = 12000000.0
    t.compute(e_o)
    print("Ephem pressure 1200mBar AzAlt: %.12f, %.12f  |  RADec: %.12f, %.12f" % (t.az, t.alt, t.ra, t.dec))
    e_o.compute_pressure()
    print("ephem_observer pressure computed: %.12f" % e_o.pressure)
    t.compute(e_o)
    print("Ephem pressure computed AzAlt: %.12f, %.12f  |  RADec: %.12f, %.12f" % (t.az, t.alt, t.ra, t.dec))
    print("--")
    back_ra, back_dec = e_o.radec_of(t.az, t.alt)
    print("Ephem AzAlt > RADec norm pressure: %s, %s  |  %.12f, %.12f" % (back_ra, back_dec, back_ra, back_dec))
    e_o.pressure = 0
    t.compute(e_o)
    back_ra, back_dec = e_o.radec_of(t.az, t.alt)
    print("Ephem AzAlt > RADec ZERO pressure: %s, %s  |  %.12f, %.12f" % (back_ra, back_dec, back_ra, back_dec))
    
    
    

if __name__ == "__main__":
    give_it_a_spin()
    #compare_pyephem_with_my_code()
    #unittest.main()
    
