#!/usr/bin/env python
# 
# Tracie Perez
#
# Miscellaneous functions needed for astrodynamics problems:
#  1) TwoBody2DEOM(state, time)
#  2) inertial2orbital(r,v,mu)
#  3) Gibbs(r1, r2, r3, mu)
#
# Reference text: Orbital Mechanics for Engineering Students
#                 by Howard D. Curtis
#                 2nd Edition
#

import sys
from constants import *
from numpy import *
from numpy import linalg as la

def TwoBody2DEOM(state, time):
  x          = state[0]
  y          = state[1]
  x_dot      = state[2]
  y_dot      = state[3]
  r_mag      = (x**2 + y**2)**0.5
  r          = array([x, y])
  r_ddot     = (-mu_E/r_mag**3)*r
  x_ddot     = r_ddot[0]
  y_ddot     = r_ddot[1]
  d_state    = zeros((4))
  d_state[0] = x_dot
  d_state[1] = y_dot
  d_state[2] = x_ddot
  d_state[3] = y_ddot
  return d_state

def inertial2orbital(r,v,mu):
  """
  This function accepts inertial position and velocity vectors and
  returns the classical orbital elements.

  This function implements Curtis' Algorithm 4.2.

  Assumption: Mass of one body is negligible compared to mass of the
  other body.

  Inputs
  r     : Inertial position vector r(t0), 3 dimensions
  v     : Inertial velocity vector v(t0), 3 dimensions
  mu    : gravitational parameter
  
  Outputs
  h_mag : magnitude of the specific angular momentum
  a     : Semimajor axis
  e_mag : Eccentricity (unitless)
  inc   : Inclination in radians
  RAAN  : Right ascension of the ascending node in radians
  omega : Argument of perigee in radians
  M0    : Mean anomaly in radians
  theta : True anomaly in radians

  """

  # Floating point precision (2**Machine_epsilon)
  epsilon = finfo(float).eps

  # Magnitude of position vector
  r_mag = la.norm(r)

  # Magnitude of velocity vector
  v_mag = la.norm(v)

  # Radial velocity
  v_radial = dot(r,v)/r_mag

  # Specific angular momentum
  h = cross(r,v)

  # Magnitude of specific angular momentum
  h_mag = la.norm(h)

  # Inclination
  inc = arccos(h[2]/h_mag)

  # Node line
  N = cross(array([0,0,1]),h)

  # Magnitude of node line
  N_mag = la.norm(N)

  # Right ascension of the ascending node
  if abs(N_mag) > epsilon:
    RAAN = arccos(N[0]/N_mag)
    if N[1] < 0:
      RAAN = 2*pi - RAAN
  else:
    RAAN = 0

  # Eccentricity
  e = (1.0/mu)*((v_mag**2-mu/r_mag)*r - r_mag*v_radial*v)

  # Magnitude of eccentricity vector
  e_mag = la.norm(e)

  # Semimajor axis
  r_p = (h_mag**2/mu)*(1/(1+e_mag))
  r_a = (h_mag**2/mu)*(1/(1-e_mag))
  a = 0.5*(r_p + r_a)

  # Argument of perigee
  if abs(N_mag) > epsilon:
    if e_mag > epsilon:
      omega = arccos(dot(N,e)/(N_mag*e_mag))
      if e[2] < 0:
        omega = 2*pi - omega
    else:
      omega = 0    
  else:
    omega = 0

  # True anomaly
  if e_mag > epsilon:
    theta = arccos(dot(e/e_mag,r/r_mag))
    if v_radial < 0:
      theta = 2*pi - theta
  else:
    cp = cross(N, r)
    if cp[2] >= 0:
      theta = arccos(dot(N/N_mag,r/r_mag))
    else:
      theta = 2*pi - arccos(dot(N/N_mag,r/r_mag))
    
  # Finding mean anomaly depends on orbit type
  if e_mag < epsilon:               # Circle
    M = theta
  elif e_mag < 1:                   # Ellipse
    # Curtis' Eqn 3.13b
    temp = ((1-e_mag)/(1+e_mag)**0.5)*tan(theta/2)
    E = 2*arctan(temp)
    if E < 0:
      E += 2*pi
    # Curtis' Eqn 3.14, Kepler's Equation
    M = E - e_mag*sin(E)
  elif e_mag == 1:                  # Parabola
    # Curtis' Eqn 3.30
    M = 0.5*tan(theta/2) + (1/6)*(tan(theta/2))**3
  else:                             # Hyperbola
    # Curtis' Eqn 3.38, F in terms of true anomaly
    temp1 = (e_mag**2-1)**0.5 * sin(theta)
    temp = temp1/(1+e_mag*cos(theta))
    # Eccentric anomaly
    F = arcsinh(temp)
    # Mean anomaly for the hyperbola, Eqn. 3.40
    M = e_mag*sinh(F)-F   

  return h_mag, a, e_mag, inc, RAAN, omega, M, theta

def Gibbs(r1, r2, r3, mu):
  """
  This function applies the Gibbs method to solve for the orbital
  elements given three geocentric position vectors r1, r2, r3. This
  algorithm is given by Curtis in Section 5.2.

  The provided position vectors much be 3D.

  The orbital elements returned are, in this order:
    a     - semimajor axis
    ecc   - eccentricity
    h     - specific angular momentum
    inc   - inclination
    RAAN  - right ascension of the ascending node
    omega - argument of perigee
    theta - true anomaly of r2
    v2    - not an orbital element, but velocity for r2
  """

  # Calculate magnitude of the position vectors
  r1_mag   = la.norm(r1)
  r2_mag   = la.norm(r2)
  r3_mag   = la.norm(r3)

  # Calculate C12 = r1 x r2, C23 = r2 x r3, C31 = r3 x r1
  C12      = cross(r1,r2)
  C23      = cross(r2,r3)
  C31      = cross(r3,r1)

  # Verify that u_r1_hat*C23_hat = 0
  u_r1_hat = r1/r1_mag
  C23_mag  = la.norm(C23)

  C23_hat  = C23/C23_mag
  test     = dot(u_r1_hat,C23_hat)
  if abs(test) > 1e-4:
    sys.stderr.write("Error on step 3 in Gibbs method function.\n")
    sys.exit(1)

  # Calculate N, D, and S using Equations 5.13, 5.14, and 5.21
  N = r1_mag*cross(r2,r3) + r2_mag*cross(r3,r1) + r3_mag*cross(r1,r2)
  D = cross(r1,r2) + cross(r2,r3) + cross(r3,r1)
  S = r1*(r2_mag-r3_mag) + r2*(r3_mag-r1_mag) + r3*(r1_mag-r2_mag)

  # Calculate v2 using Equation 5.22
  N_mag = la.norm(N)
  D_mag = la.norm(D)
  denom = N_mag*D_mag
  v2    = sqrt(mu/denom)*((1/r2_mag)*cross(D,r2) + S)

  # Use r2 and v2 to compute the orbital elements via Algorithm 4.2
  [h,a,ecc,inc,RAAN,omega,_,theta] = inertial2orbital(r2,v2,mu_E)

  return a, ecc, h, inc, RAAN, omega, theta, v2