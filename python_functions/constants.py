#!/usr/bin/env python
# 
# Tracie Perez
#
# Constants for use in the orbit simulations
#
# These are purposefully not integers to avoid rounding errors when
# using them in equations with floats. To have Python do proper
# division, add this line to top of programs:
# from __future__ import division\
#

mu_E_m  = 3.986e14  # m^3/s^2, standard gravitational param. of Earth
mu_E    = 398600.0  # km^3/s^2, ditto
r_E     = 6378.0    # km, radius of Earth
