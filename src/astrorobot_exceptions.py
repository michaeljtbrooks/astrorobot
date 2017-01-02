#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    Exceptions
"""

from __future__ import unicode_literals

class InitError(Exception):
    """
    A class to represent incorrect setting up of your situation 
    """
    pass

class DownloadingError(Exception):
    """
    When a download / http fetch fails.
    """
    pass