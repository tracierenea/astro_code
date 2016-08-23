#!/usr/bin/env python
# 
# Tracie Perez
#
# Plot the synthetic measurements and state estimate results.
#

from numpy import *
import matplotlib.pyplot as plt
import os
import sys

### Check that one argument was passed
num_args = len(sys.argv)
if num_args != 2:
  print "Error in create_plots_MC.py:"
  exit('  --> Provide argument of the test case ran: 1, 2, or 3')
test_case = int(sys.argv[1])

#### Read in data from file
meas_data = loadtxt('measurement_data.txt',
                    skiprows=1)              # time, y_true, y_meas
states    = loadtxt('sat_states.txt')        # size m x 17. see [1]

with open('measurement_data.txt', 'r') as f:
    first_line = f.readline().rstrip()
x0_guess = array(map(float, first_line.split()))

#### General parameters/settings
r_Titan = 2575.5  # km, radius of Titan

#### Needed to put image of Titan on figure
# The full path to the image file must be specified for the function
# get_sample_data, otherwise it assumes the path is relative to the
# mpl-data/sample_data directory
Titan_image_file = os.getcwd() + "/Titan_cropped.png"
# Read in image and plot Titan
image = plt.imread(Titan_image_file)

#################### Plot 1: measurement data #####################
fig = plt.gcf()
plt.plot(meas_data[:,0]/60.0, meas_data[:,2], 'o', 
         label='Measurements', alpha=0.6, markerfacecolor='blue',
         markeredgecolor='blue')
plt.plot(meas_data[:,0]/60.0, meas_data[:,1], 'k-', label='Truth')
plt.xlabel('Time (minutes)',     fontsize=16)
plt.ylabel('Doppler Shift (Hz)', fontsize=16)
plt.title('Measurements Created for Simulation', fontsize=16)
plt.legend(shadow=True, fontsize=16, loc='best')
plt.savefig("figure1_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 2: satellite positions #################
fig = plt.gcf()
# Add image of Titan
plt.imshow(image, extent=[-r_Titan,r_Titan,-r_Titan,r_Titan])
# TODO: Would be nice to add the initial guess to the figure, but
#       that isn't contained in either file
plt.plot(states[:,1], states[:,2], 'ko', label='Mothersat')
plt.plot(states[:,5], states[:,6], 's', label='Femtosat (truth)',
         markerfacecolor='green', markeredgecolor='green')
plt.plot(states[:,9], states[:,10], 's',
         label='Femtosat (NLS+EKF estimate)', markerfacecolor='red',
         markeredgecolor='red')
plt.legend(shadow=True, fontsize=14, loc='upper left')
# Start and end points
end = states.shape[0] - 1
plt.text(states[0,1]-1500,   states[0,2]-1500,  "Start", fontsize=14)
plt.text(states[end,1]-1500, states[end,2]+500, "End",   fontsize=14)
plt.text(states[0,9]+500,    states[0,10],      "Start", fontsize=14)
plt.text(states[end,9]+500,  states[end,10],    "End",   fontsize=14)
plt.axis('equal')
plt.xlabel('x (km)', fontsize=16)
plt.xlabel('y (km)', fontsize=16)
plt.title('Estimated and True Femtosatellite Position', fontsize=16)
plt.savefig("figure2_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 3: measurement availability ############
# Add image of Titan
plt.imshow(image, extent=[-r_Titan,r_Titan,-r_Titan,r_Titan])

n = states.shape[0]
time_meas = meas_data[:,0]

# Only plot a point every 30 seconds
for index in arange(0,n,30):
  plt.plot(states[:,1], states[:,2], 'ko')
  # Check if this time is in the measurement time column. If so, plot
  # in green; otherwise, plot in red.
  temp_index = where(time_meas == index)[0]
  if temp_index.size:
    plt.plot(states[index,5], states[index,6], 's', 
             markerfacecolor='green', markeredgecolor='green')
  else:
    plt.plot(states[index,5], states[index,6], 's',
             markerfacecolor='red', markeredgecolor='red')
plt.text(states[0,1],       states[0,2]-350,    "Start", fontsize=14)
plt.text(states[end,1]-100, states[end,2]-400,  "End",   fontsize=14)
plt.text(states[0,9]+100,   states[0,10]-400,   "Start", fontsize=14)
plt.text(states[end,9]+100, states[end,10]-400, "End",   fontsize=14)
plt.axis('equal')
plt.xlabel('x (km)', fontsize=16)
plt.xlabel('y (km)', fontsize=16)
plt.title('Measurement Availability\nGreen = possible, Red = not possible due to signal blockage', fontsize=18)
plt.savefig("figure3_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plots 4 & 5: femtosat position error ###########
# Residual  = estimate - truth
x_resids_KF = states[:,9]  - states[:,5]
y_resids_KF = states[:,10] - states[:,6]
resids_KF   = (x_resids_KF**2 + y_resids_KF**2)**0.5
plt.plot(states[:,0]/60.0, resids_KF, 'k.')
plt.grid('on')
plt.xlabel('Time (minutes)', fontsize=16)
plt.ylabel('Error (km)',     fontsize=16)
plt.title('Femtosatellite Position Estimate Error', fontsize=16)
plt.savefig("figure4_case" + str(test_case) + ".png")
plt.show()
fig.clear()

# Now we also want to see what the error would have been had we just
# propagated the initial guess, without filtering any measurements.
# So start by re-plotting the same KF estimate residuals.
plt.plot(states[:,0]/60.0, resids_KF, 'k-', label='KF Residuals')
plt.grid('on')
plt.xlabel('Time (minutes)', fontsize=16)
plt.ylabel('Error (km)',     fontsize=16)
plt.title('Femtosatellite Position Estimate Error', fontsize=16)
# Next, propagate the initial guess
def Two_body_EOM(state, time):
  G          = 6.6742e-11 / 1000**3;# km^3/kg*s^2, Univ grav constant
  m_Titan    = 1.3452e23;           # kg, mass of Titan
  mu         = G*m_Titan;           # km^3/s^2, Titan's grav. param
  x          = state[0]
  y          = state[1]
  x_dot      = state[2]
  y_dot      = state[3]

  r_mag      = (x**2 + y**2)**0.5
  r          = array([x, y])
  r_ddot     = (-mu/r_mag**3)*r
  x_ddot     = r_ddot[0]
  y_ddot     = r_ddot[1]
  d_state    = zeros((4))
  d_state[0] = x_dot
  d_state[1] = y_dot
  d_state[2] = x_ddot
  d_state[3] = y_ddot
  return d_state
from scipy.integrate import odeint
time_array     = states[:,0]
results        = odeint(Two_body_EOM, x0_guess, time_array)
x_2body        = results[:,0]
y_2body        = results[:,1]
x_resids_2body = x_2body - states[:,5]
y_resids_2body = y_2body - states[:,6]
resids_2body   = (x_resids_2body**2 + y_resids_2body**2)**0.5
initial_error  = ((x0_guess[0]-states[0,5])**2 + 
                  (x0_guess[1]-states[0,6])**2)**0.5
print "True initial position : %.3f" % states[0,5], " %.3f" % \
      states[0,6]
print "Initial position guess: %.3f" % x0_guess[0], " %.3f" % \
      x0_guess[1]
print "Initial guess error   : %.3f" % initial_error
plt.plot(states[:,0]/60.0, resids_2body, 'r-', 
         label='Guess Propagation Residuals')
plt.legend(shadow=True, fontsize=14, loc='best')
plt.savefig("figure5_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 6: femtosat velocity error ###########
# Residual  = estimate - truth
x_dot_resids_KF = states[:,11] - states[:,7]
y_dot_resids_KF = states[:,12] - states[:,8]
vel_resids_KF   = (x_dot_resids_KF**2 + y_dot_resids_KF**2)**0.5
plt.plot(states[:,0]/60.0, vel_resids_KF, 'k.')
plt.grid('on')
plt.xlabel('Time (minutes)', fontsize=16)
plt.ylabel('Error (km/sec)',     fontsize=16)
plt.title('Femtosatellite Velocity Estimate Error', fontsize=18)
plt.savefig("figure6_case" + str(test_case) + ".png")
plt.show()
fig.clear()



#################### Plot 7: state error & 3sig bounds ###########
# Error  = estimate - truth
time_array = states[:,0]/60.0
fig = plt.figure()
ax1 = fig.add_subplot(211)
x_KF_error = states[:,9]  - states[:,5]
ax1.plot(time_array, x_KF_error, 'k-', label="x error")
ax1.plot(time_array,  states[:,13], 'g--', label =  "P_xx forward")
ax1.plot(time_array, -states[:,13], 'g--', label = "-P_xx forward")
ax1.plot(time_array,  states[:,15], 'b--', label =  "P_xx backward")
ax1.plot(time_array, -states[:,15], 'b--', label = "-P_xx backward")
plt.legend(shadow=True, fontsize=14, loc='best')
ax2 = fig.add_subplot(212)
y_KF_error = states[:,10]  - states[:,6]
ax2.plot(time_array, y_KF_error, 'k-', label="y error")
ax2.plot(time_array,  states[:,14], 'g--', label =  "P_yy forward")
ax2.plot(time_array, -states[:,14], 'g--', label = "-P_yy forward")
ax2.plot(time_array,  states[:,16], 'b--', label =  "P_yy backward")
ax2.plot(time_array, -states[:,16], 'b--', label = "-P_yy backward")
plt.legend(shadow=True, fontsize=14, loc='best')
plt.savefig("figure7_case" + str(test_case) + ".png")
plt.show()
fig.clear()
####################################################################

# footnote [1]
#   The columns of the array states will be, from left to right:
#   1) time of the state           index: 0
#   2) r_x     of mothersat (truth)       1
#   3) r_y     of mothersat (truth)       2
#   4) r_x_dot of mothersat (truth)       3
#   5) r_y_dot of mothersat (truth)       4
#   6) r_x     of femsat    (truth)       5
#   7) r_y     of femsat    (truth)       6
#   8) r_x_dot of femsat    (truth)       7
#   9) r_y_dot of femsat    (truth)       8
#  10) r_x     of femsat    (estimate)    9
#  11) r_y     of femsat    (estimate)    10
#  12) r_x_dot of femsat    (estimate)    11
#  13) r_y_dot of femsat    (estimate)    12
#  14) P_xx covariance, forward pass      13
#  15) P_yy covariance, forward pass      14
#  16) P_xx covariance, backwards pass    15
#  17) P_yy covariance, backwards pass    16
#
#
# footnote [2]
#   To plot a blue circle instead of Earth image, use:
#     Earth_outline = plt.Circle((0,0), r_Titan, color='blue')
#     fig.gca().add_artist(Earth_outline)