#!/usr/bin/env python
# 
# Tracie Perez
#
# A QA file to check that the functions in astro_functions are
# working correctly.
#
# Examples used are from: Orbital Mechanics for Engineering Students
#                         by Howard D. Curtis
#                         2nd Edition
#

from constants import *
from astro_functions import *
from numpy import *
from numpy import linalg as la

pass_small  = 0.1   # Degrees
pass_bigger = 4;    # units of h, a

###########################
# Test inertial2orbital function by using values from Example 4.3 on
# page 212.
print "Test #1) Example 4.3, results from inertial2orbital function:"
r         = array([ -6045, -3490,  2500])
v         = array([-3.457, 6.618, 2.533])

[h_mag,a,e_mag,inc,RAAN,omega,M,theta] = inertial2orbital(r, v, mu_E)

inc_deg   = inc*180/pi
RAAN_deg  = RAAN*180/pi
omega_deg = omega*180/pi
theta_deg = theta*180/pi

print "h     - specific angular momentum:   %.1f" % h_mag,
if abs(h_mag - 58310) < pass_bigger:
  print "\tPASSED, delta = %.1f" % abs(h_mag - 58310)
else:
  print "\tFAILED, truth = 58,310 km^2/s"
print "a     - semimajor axis:              %.1f" % a,
if abs(a - 8788) < pass_bigger:
  print "\tPASSED, delta = %.1f" % abs(a - 8788)
else:
  print "\tFAILED, truth = 8788 km"
print "ecc   - eccentricity:                %.3f" % e_mag,
if abs(e_mag - 0.1712) < pass_small:
  print "\tPASSED, delta = %.3f" % abs(e_mag - 0.1712)
else:
  print "\tFAILED, truth = 0.1712"
print "inc   - inclination:                 %.1f" % inc_deg,
if abs(inc_deg - 153.2) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(inc_deg - 153.2)
else:
  print "\tFAILED, truth = 153.2 deg"
print "RAAN  - right asc. of the asc. node: %.1f" % RAAN_deg,
if abs(RAAN_deg - 255.3) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(RAAN_deg - 255.3)
else:
  print "\tFAILED, truth = 255.3 deg"
print "omega - argument of perigee:         %.3f" % omega_deg,
if abs(omega_deg - 20.07) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(omega_deg - 20.07)
else:
  print "\tFAILED, truth = 20.07 deg"
print "theta - true anomaly:                %.1f" % theta_deg,
if abs(theta_deg - 28.45) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(theta_deg - 28.45)
else:
  print "\tFAILED, truth = 28.45 deg"
print

###########################
# Test inertial2orbital function by using values from Problem 4.3 on
# page 250.
print "Test #2) Problem 4.3, results from inertial2orbital function:"
r         = array([  2615,   15881,  3980])
v         = array([-2.767, -0.7905, 4.980])

[h_mag,a,e_mag,inc,RAAN,omega,M,theta] = inertial2orbital(r, v, mu_E)

inc_deg   = inc*180/pi
RAAN_deg  = RAAN*180/pi
omega_deg = omega*180/pi
theta_deg = theta*180/pi

print "h     - specific angular momentum:   %.1f" % h_mag,
if abs(h_mag - 95360) < pass_bigger:
  print "\tPASSED, delta = %.1f" % abs(h_mag - 95360)
else:
  print "\tFAILED, truth = 95,360 km^2/s"
print "ecc   - eccentricity:                %.3f" % e_mag,
if abs(e_mag - 0.3760) < pass_small:
  print "\tPASSED, delta = %.3f" % abs(e_mag - 0.3760)
else:
  print "\tFAILED, truth = 0.3760"
print "inc   - inclination:                 %.1f" % inc_deg,
if abs(inc_deg - 63.95) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(inc_deg - 63.95)
else:
  print "\tFAILED, truth = 63.95 deg"
print "RAAN  - right asc. of the asc. node: %.1f" % RAAN_deg,
if abs(RAAN_deg - 73.71) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(RAAN_deg - 73.71)
else:
  print "\tFAILED, truth = 73.71 deg"
print "omega - argument of perigee:         %.3f" % omega_deg,
if abs(omega_deg - 15.43) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(omega_deg - 15.43)
else:
  print "\tFAILED, truth = 15.43 deg"
print "theta - true anomaly:                %.1f" % theta_deg,
if abs(theta_deg - 0.06764) < pass_small:
  print "\tPASSED, delta = %.3f deg\n" % abs(theta_deg - 0.06764)
else:
  print "\tFAILED, truth = 0.06764 deg\n"

###########################
# Test inertial2orbital function by using values from Problem 4.4 on
# page 250.
print "Test #3) Problem 4.4, results from inertial2orbital function:"
r         = array([0,      0,   12670])
v         = array([0, -3.874, -0.7905])

[h_mag,a,e_mag,inc,RAAN,omega,M,theta] = inertial2orbital(r, v, mu_E)

inc_deg   = inc*180/pi
RAAN_deg  = RAAN*180/pi
omega_deg = omega*180/pi
theta_deg = theta*180/pi

print "h     - specific angular momentum:   %.1f" % h_mag,
if abs(h_mag - 49080) < pass_bigger:
  print "\tPASSED, delta = %.1f" % abs(h_mag - 49080)
else:
  print "\tFAILED, truth = 49,080 km^2/s"
print "ecc   - eccentricity:                %.3f" % e_mag,
if abs(e_mag - 0.5319) < pass_small:
  print "\tPASSED, delta = %.3f" % abs(e_mag - 0.5319)
else:
  print "\tFAILED, truth = 0.5319"
print "inc   - inclination:                 %.1f" % inc_deg,
if abs(inc_deg - 90) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(inc_deg - 90)
else:
  print "\tFAILED, truth = 90 deg"
print "RAAN  - right asc. of the asc. node: %.1f" % RAAN_deg,
if abs(RAAN_deg - 90) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(RAAN_deg - 90)
else:
  print "\tFAILED, truth = 90 deg"
print "omega - argument of perigee:         %.3f" % omega_deg,
if abs(omega_deg - 259.50) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(omega_deg - 259.50)
else:
  print "\tFAILED, truth = 259.50 deg"
print "theta - true anomaly:                %.1f" % theta_deg,
if abs(theta_deg - 190.50) < pass_small:
  print "\tPASSED, delta = %.3f deg\n" % abs(theta_deg - 190.50)
else:
  print "\tFAILED, truth = 190.50 deg\n"

###########################
# Test inertial2orbital function by using values from Problem 4.5 on
# page 250.
print "Test #4) Problem 4.5, results from inertial2orbital function:"
r         = array([6472.7, -7470.8, -2469.8])
v         = array([3.9914,  2.7916, -3.2948])

[h_mag,a,e_mag,inc,RAAN,omega,M,theta] = inertial2orbital(r, v, mu_E)

inc_deg   = inc*180/pi
RAAN_deg  = RAAN*180/pi
omega_deg = omega*180/pi
theta_deg = theta*180/pi

print "h     - specific angular momentum:   %.1f" % h_mag,
if abs(h_mag - 58461) < pass_bigger:
  print "\tPASSED, delta = %.1f" % abs(h_mag - 58461)
else:
  print "\tFAILED, truth = 58,461 km^2/s"
print "ecc   - eccentricity:                %.3f" % e_mag,
if abs(e_mag - 0.2465) < pass_small:
  print "\tPASSED, delta = %.3f" % abs(e_mag - 0.2465)
else:
  print "\tFAILED, truth = 0.2465"
print "inc   - inclination:                 %.1f" % inc_deg,
if abs(inc_deg - 35) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(inc_deg - 35)
else:
  print "\tFAILED, truth = 35 deg"
print "RAAN  - right asc. of the asc. node: %.1f" % RAAN_deg,
if abs(RAAN_deg - 110) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(RAAN_deg - 110)
else:
  print "\tFAILED, truth = 110 deg"
print "omega - argument of perigee:         %.3f" % omega_deg,
if abs(omega_deg - 75) < pass_small:
  print "\tPASSED, delta = %.3f deg" % abs(omega_deg - 75)
else:
  print "\tFAILED, truth = 75 deg"
print "theta - true anomaly:                %.1f" % theta_deg,
if abs(theta_deg - 130) < pass_small:
  print "\tPASSED, delta = %.3f deg\n" % abs(theta_deg - 130)
else:
  print "\tFAILED, truth = 130 deg\n"

###########################
# Test Gibbs function by using values from Example 5.1 on page 261.
print "Test #5) Example 5.1, results from Gibbs function:"
r1        = array([-294.32, 4265.1, 5986.7])
r2        = array([-1365.5, 3637.6, 6346.8])
r3        = array([-2940.3, 2473.7, 6555.8])

[a, ecc, h, inc, RAAN, omega, theta, v2] = Gibbs(r1, r2, r3, mu_E)

inc_deg   = inc*180/pi
RAAN_deg  = RAAN*180/pi
omega_deg = omega*180/pi
theta_deg = theta*180/pi

print "a     - semimajor axis:              %.1f" % a,
if abs(a - 8000) < pass_bigger:
  print "\tPASSED"
else:
  print "\tFAILED"
print "ecc   - eccentricity:                %.3f" % ecc,
if abs(ecc - 0.1) < pass_small:
  print "\tPASSED"
else:
  print "\tFAILED"
print "h     - specific angular momentum:   %.3f" % h
print "inc   - inclination:                 %.1f" % inc_deg,
if abs(inc_deg - 60) < pass_small:
  print "\tPASSED"
else:
  print "\tFAILED"
print "RAAN  - right asc. of the asc. node: %.1f" % RAAN_deg,
if abs(RAAN_deg - 40) < pass_small:
  print "\tPASSED"
else:
  print "\tFAILED"
print "omega - argument of perigee:         %.3f" % omega_deg,
if abs(omega_deg - 30) < pass_small:
  print "\tPASSED"
else:
  print "\tFAILED"
print "theta - true anomaly for r2:         %.1f" % theta_deg,
if abs(theta_deg - 50) < pass_small:
  print "\tPASSED\n"
else:
  print "\tFAILED\n"

###########################
# Test Gibbs function by using values from Problem 5.1 on page 312.
print "Test #6) Problem 5.1, results from Gibbs function:"
r1        = array([5887, -3520, -1204])
r2        = array([5572, -3457, -2376])
r3        = array([5088, -3289, -3480])

[a, ecc, h, inc, RAAN, omega, theta, v2] = Gibbs(r1, r2, r3, mu_E)

inc_deg   = inc*180/pi
RAAN_deg  = RAAN*180/pi
omega_deg = omega*180/pi
theta_deg = theta*180/pi
print "v2  : %.3f km/sec" % la.norm(v2),
if abs(la.norm(v2) - 7.59):
  print "\t\t\t\tPASSED"
else:
  print "\t\t\t\tFAILED"

r_p = (1/mu_E)*(h**2)*(1/(1+ecc)) # Eqn. 2.50 on page 81
z_p = r_p - r_E
print "z_p : %.3f km" % z_p,
if abs(z_p - 567):
  print "\t\t\t\tPASSED\n"
else:
  print "\t\t\t\tFAILED\n"

###########################
# Test Gibbs function. Use the elliptical orbit given in Example 2.7
print "Test #7) Example 2.7 elliptical orbit, results from Gibbs function:"
z_p = 400.0                       # km
z_a = 4000.0                      # km
r_a = z_a + r_E                   # km
r_p = z_p + r_E                   # km
e   = (r_a - r_p)/(r_a + r_p)     # eccentricity
a   = 0.5*(r_p + r_a)             # semimajor axis
b   = a*sqrt(1-e**2)              # semiminor axis
r1  = array([ r_p, 0, 0])         # position vector at periapsis
r2  = array([-r_a, 0, 0])         # position vector at apoapsis
r3  = array([-a*e, b, 0])         # look at Fig. 2.18 on pg 89

[a, ecc, h, inc, RAAN, omega, theta, v2] = Gibbs(r1, r2, r3, mu_E)

theta_deg = theta*180/pi

print "a     - semimajor axis:              %.1f" % a,
if abs(a - 8578) < pass_bigger:
  print "\tPASSED"
else:
  print "\tFAILED"
print "ecc   - eccentricity:                %.3f" % ecc,
if abs(ecc - 0.2098) < pass_small:
  print "\tPASSED"
else:
  print "\tFAILED"
print "h     - specific angular momentum:   %.3f" % h,
if abs(h - 57172) < pass_bigger:
  print "\tPASSED, delta = %.3f" % abs(h - 57172)
else:
  print "\tFAILED, truth = 57,172 km^2/s"
print "theta - true anomaly for r2:         %.1f" % theta_deg,
if abs(theta_deg - 180) < pass_small:
  print "\tPASSED"
else:
  print "\tFAILED"
print "v at apoapsis: %.3f km/sec" % la.norm(v2),
if abs(la.norm(v2) - 5.509) < pass_small:
  print "\t\t\tPASSED"
else:
  print "\t\t\tFAILED"
