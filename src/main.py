#!/usr/bin/python
# -*- coding: utf8 -*-

"""
AstroRobot
    
    Controls all aspects of an astrophotography scope and camera for obtaining images.
    
    Elements to be controlled:
    
        Initially:
            - Nikon Camera Shutter
            - Scope position (tapping stepper motor controller)
                * uses "ballistic" method of: measure motor step commands and assuming each one results in a true step
        
        Later:
            - Filter wheel
            - Image acquisition (inc flipping back if Meridian flip has occurred) 
            - Cooled CMOS camera capture
            - Autoguiding
            - Scope focus
    
    Nomenclature:
        - azalt / altaz = azimuth & altitude co-ordinates
        - equatorial_mount / mount = mount co-ordinates = the Earth's Lat and Lon the scope would be perpendicular to if it were mounted at the centre of the Earth
        - celestial = Right Ascension & Declination co-ordinates (will drift on non-tracking!)
    
    
    Dependencies / Other bits to add in:
        - PyNOVAS = Astronomical functions and applications (http://pynovas.sourceforge.net/)
        
    
    
    @author:   Dr Michael J T Brooks
    @date:     2016-02-24
    @version:  20160224  
"""


"""
    Components needed:
    
        main.py = the main hardware routine (i.e. represents the Pi). Needs to be an infallable loop!
            ? Do we use threading to partition the different jobs, or shall we keep it simple do it all in one thread?
                If we're having multiple listeners, this suggests we should be THREADING
        
        
        
        telescope.py = the interface to the telescope, including controller
            Needs to:
                - cope with the "Meridian flip"
                - have a concept of "forbidden zones" which is positions the scope CANNOT enter (e.g. positions which will smash the scope into the floor)
                - knows that the current sun position is in the forbidden zone UNLESS the scope is a solar scope
                
            
            class EQ5
                Properties:
                    location                     Where on Earth the scope is positioned! - Could pull this in from GPS data?
                    mount                        Defines the type of mount (AzAlt / equatorial)... in future could use to model completely new axes!
                    scope_side                   W or E: Which side of the mount the scope is. When looking at an object crossing the meridian to the West, needs to be on the East to avoid pier clonkage
                    current_position             Where it is pointing (a Position object) - only motor movements can change this!
                    target                       the current target (returns SkyObject or None)
                    on_target                    True / False - If the scope thinks it is looking at its target!
                    tracking                     True / False - if the scope is tracking
                    solar                        True / False - is this a solar scope? (if so, will track the sun!)
                    initialised_time             Precise UTC time the telescope was initialised (needs to happen after NTP fix!)
                      
                
                Methods:
                    position_at(<UTC time>)      Returns the position the scope will be at the given time. If tracking, RA/Dec will stay constant, otherwise AzAlt + LatLon will remain constant
                    move_time(PositionDelta)     Works out how long it will take the motors to move through the given range, INCLUDING accounting for Meridian flip
                    _move_time(PositionDelta)    (Sub-method which calculates move time EXCLUDING the meridian flip!)
                    go_to(SkyObject | Position)  Makes the scope move to the given target
                                                 Calculates the move time, and adds the relevant amount of Earth rotation onto it so the motors compensate for any
                                                 earth rotation during their movement.
                    
                    motor_listener()             Listens for motor movements and updates the current_position according to the motor movements
                                                 it detects. Thus steady tracking will be detected as an increasing RA, BUT this will get corrected
                                                 using the rotation_adjust method
                    rotation_adjust()            Calculates how much Earth rotation has occured since the last check and subtracts away the
                                                 relevant angle in RA from the telescope's current position. The motor_listener will have
                                                 nudged the scope on in RA if in tracking mode, thus this will result in a zero change if
                                                 tracking correctly
                                 
                    
        
        
        
        celestialsphere.py = represents the sky and where objects are, given the time of day and the location, talks to the relevant planetarium database
        
            
            class BasePosition = a clever object that can map between Az-alt, RA Dec and our own scope polar py
                Properties:
                    azalt               Azimuth & Altitude co-ordinates
                    equatorial_mount    Lat Lon perpendicular mount co-ordinates
                    celestial           Right Ascension & Declination sky co-ordinates (moves relative to earth)
                    true_at_time        UTC timestamp that this Position applies to
                
                Methods:
                    get_position_at(<time>)    Returns a new Position for the specified time 
                    calc_subtract()     if second object is a PositionDelta, will give the new Position
                        _calc_sub_pos() if second object is a Position, will give a PositionDelta to indicate moving from one position to the other
                        _calc_sub_posd()
                    calc_add()          if second object is a PositionDelta, will give the new Position
                        _calc_add_pos() if second object is a Position, this makes no sense so will barf error
        
            class PositionTerrestrial  = A position where AzAlt / LatperpLonperp is fixed and RA,Dec varies with time, applies to scope
                Properties
                    tracked = False
            
            class PositionCelestial  = A position where RA,Dec is fixed and AzAlt / LatperpLonperp is varies with time, applies to sky objects 
                Properties
                    tracked = False
        
            
            class PositionDelta = an analogue to Position which represents the change in the three grid py
                Properties
                    azalt               Change in Azimuth & Altidude co-ordinates
                    equatorial_mount    Change in Lat Lon perpendicular mount co-ordinates
                    celestial           Change in Right Ascension & Declination sky co-ordinates (moves relative to earth)
                    crosses_meridian    True / False / None - If the PositionDelta has been calculated from two Position objects, and they straddle the meridian
                                        then this will be True, indicating that might not be safe to move directly from one to the other. None = unknown 
                
                Methods
                    move_time_celestial    How long it will take motors to move through this PositionDelta, accounting for earth rotation
                    move_time_terrestrial  How long it will take motors to move through this PositionDelta, ignoring earth rotation (e.g. to hit AzAlt targets)
                    _move_time_uncorrected sub method used in both. 
        
            
            class SkyObject = an object in the sky. Given your current GPS location and time of day, will convert between RA,Dec to Alt,Az 
                Properties:
                    name               the object's colloquial name (if it has one)
                    messier            the object's messier name (if it has one)
                    ngc                the object's NGC name (if it has one)
                    py        the location on the celestial sphere in RA,Dec
                    position           a Position instance showing where the object resides in all three axes
                    celestial_type     Whether a star / galaxy / nebula / planet / 
                    apparent_magnitude          
                    
                Methods:    
                    is_above_horizon()    Whether this object is above your 0deg horizon
                    visible()             Whether this object is actually visible given magnitude, current moon phase, 
                
        
        cameras.py = the interfaces to various cameras
            class DSLRshutter = simple shutter controls for DSLRs e.g. Nikon
            class CoolCMOS = interface to whole cooled CMOS camera 
        
"""