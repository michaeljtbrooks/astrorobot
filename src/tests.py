#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot

    Tests
"""
import unittest
#
from coordinates import *


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


if __name__ == "__main__":
    unittest.main()