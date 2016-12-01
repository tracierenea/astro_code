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
  print "Error in create_plots.py:"
  exit('  --> Provide argument of the test case ran: 1, 2, or 3')
test_case = int(sys.argv[1])

#### Read in data from file. F denotes the forward-only KF results;
#### FB denotes the forward-backwards results.
meas_data = loadtxt('measurement_data.txt',
                    skiprows=1)              # time, y_true, y_meas
states    = loadtxt('sat_states.txt',
                    skiprows=1)              # size m x 17. see [1]
with open('measurement_data.txt', 'r') as f:
    first_line = f.readline().rstrip()
x0_guess = array(map(float, first_line.split()))

with open('sat_states.txt', 'r') as f:     #### same in both cases?
    first_line = f.readline().rstrip()
x0_fem_truth = array(map(float, first_line.split()))

#### General parameters/settings
r_Titan = 2575.5  # km, radius of Titan

#### Needed to put image of Titan on figure
# The full path to the image file must be specified for the function
# get_sample_data, otherwise it assumes the path is relative to the
# mpl-data/sample_data directory
Titan_image_file = os.getcwd() + "/Titan_cropped.png"
# Read in image and plot Titan
image = plt.imread(Titan_image_file)

# Change the default font properties
AxesFont  = {'family': 'serif', 'size': 16 }
TitleFont = {'family': 'serif', 'size': 18 }
LegFont   = {'family': 'serif', 'size': 12 }
x_label   = 'Time (minutes)'

# Save/label data for further use
x_momsat  = states[:,1]  # mothersat true x-coordinates
y_momsat  = states[:,2]  # mothersat true y-coordinates
x_truth   = states[:,5]  # femtosat true x-coordinate
y_truth   = states[:,6]  # femtosat true y-coordinate
xd_truth  = states[:,7]  # femtosat true x_dot
yd_truth  = states[:,8]  # femtosat true y_dot
x_est     = states[:,9]  # femtosat estimated x-coord, no check
y_est     = states[:,10] # femtosat estimated y-coord, no check
xd_est    = states[:,11] # femtosat estimated x_dot, no check
yd_est    = states[:,12] # femtosat estimated y_dot, no check

#################### Plot 1: measurement data #####################
y_label = 'Doppler Shift (kHz)'
title   = ' Measurements Created for Simulation'

fig = plt.gcf()
plt.plot(meas_data[:,0]/60.0, meas_data[:,2]/1000.0, 'o', 
         label='Measurements', alpha=0.6, markerfacecolor='blue',
         markeredgecolor='blue')
plt.plot(meas_data[:,0]/60.0, meas_data[:,1]/1000.0, 'k-', 
         label='Truth')
plt.xticks(fontsize=14, family='serif')
plt.yticks(fontsize=14, family='serif')
plt.xlabel(x_label, fontdict = AxesFont)
plt.ylabel(y_label, fontdict = AxesFont)
plt.title(title,    fontdict = TitleFont)
plt.legend(prop=LegFont, shadow='True', loc='upper left')
plt.savefig("figure1_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 2: satellite positions #################
title = 'Estimated and True Femtosatellite Position'

# For adding text to the plot
end = states.shape[0] - 1
if test_case == 1:
  mom_start_x = states[0,1] + 100
  mom_start_y = states[0,2] - 100 
  mom_end_x   = states[end,1] + 100
  mom_end_y   = states[end,2]
  # just make them on top of eachother - too close in this case
  fem_start_x = mom_start_x
  fem_start_y = mom_start_y
  fem_end_x   = mom_end_x
  fem_end_y   = mom_end_y
elif test_case == 2:
  mom_start_x = states[0,1] + 100
  mom_start_y = states[0,2] - 100 
  mom_end_x   = states[end,1] + 100
  mom_end_y   = states[end,2] + 100
  fem_start_x = states[0,9] + 100
  fem_start_y = states[0,10]
  fem_end_x   = states[end,9] + 100
  fem_end_y   = states[end,10]
else:
  mom_start_x = states[0,1]   + 150
  mom_start_y = states[0,2]   - 100 
  mom_end_x   = states[end,1] + 100
  mom_end_y   = states[end,2] - 100
  fem_start_x = states[0,9]   + 150
  fem_start_y = states[0,10]  - 100
  fem_end_x   = states[end,9] + 200
  fem_end_y   = states[end,10]- 100

# Add image of Titan
fig = plt.gcf()
plt.imshow(image, extent=[-r_Titan,r_Titan,-r_Titan,r_Titan])
# TODO: Would be nice to add the initial guess to the figure, but
#       that isn't contained in either file
plt.plot(x_momsat, y_momsat, 'ko', label='Mothersat')
plt.plot(x_truth,  y_truth, 's', label='Femtosat (truth)',
         markerfacecolor='green', markeredgecolor='green')
plt.plot(x_est, y_est, 's', label='Femtosat (NLS+EKF estimate)',
         markerfacecolor='red', markeredgecolor='red')
plt.legend(prop=LegFont, shadow='True', loc='upper left')
# Label which end of the trajectory is start/end
plt.text(mom_start_x, mom_start_y,  "Start", fontdict=LegFont)
plt.text(mom_end_x,   mom_end_y,    "End",   fontdict=LegFont)
plt.text(fem_start_x, fem_start_y,  "Start", fontdict=LegFont)
plt.text(fem_end_x,   fem_end_y,    "End",   fontdict=LegFont)
plt.axis('equal')
plt.xticks(fontsize=14, family='serif')
plt.yticks(fontsize=14, family='serif')
plt.xlabel('x (km)', fontdict = AxesFont)
plt.ylabel('y (km)', fontdict = AxesFont)
plt.title( title,    fontdict = TitleFont)
plt.savefig("figure2_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 3: measurement availability ############
title = 'Measurement Availability'

# Add image of Titan. See footnote [2] for how to plot circle instead
plt.imshow(image, extent=[-r_Titan,r_Titan,-r_Titan,r_Titan])

n = states.shape[0]
time_meas = meas_data[:,0]

# Only plot a point every 30 seconds
for index in arange(0,n,30): 
  # Check if this time is in the measurement time column. 
  temp_index = where(time_meas == index)[0]

  # If measurement was available, plot in green. Otherwise, plot in
  # red.
  if temp_index.size:
    # Plot the mothersat's true position at this epoch
    plt.plot(x_momsat[index], y_momsat[index], 's', 
             markerfacecolor='green', markeredgecolor='black')

    # Plot the femtosat's true position at this epoch
    plt.plot(x_truth[index], y_truth[index], 's', 
             markerfacecolor='green', markeredgecolor='green')
  else:
    # Plot the mothersat's true position at this epoch
    plt.plot(x_momsat[index], y_momsat[index], 's', 
             markerfacecolor='red', markeredgecolor='black')

    # Plot the femtosat's true position at this epoch
    plt.plot(x_truth[index], y_truth[index], 's', 
             markerfacecolor='red', markeredgecolor='red')
plt.text(mom_start_x, mom_start_y,  "Start", fontdict=LegFont)
plt.text(mom_end_x,   mom_end_y,    "End",   fontdict=LegFont)
plt.text(fem_start_x, fem_start_y,  "Start", fontdict=LegFont)
plt.text(fem_end_x,   fem_end_y,    "End",   fontdict=LegFont)

# This is just to show that for some reason, the image is being
# plotted slightly too small...
# plt.plot(-r_Titan, 0, 'k*')
# plt.plot(r_Titan, 0, 'k*')
# plt.plot(0, r_Titan, 'k*')
# plt.plot(0, -r_Titan, 'k*')

plt.xticks(fontsize=14, family='serif')
plt.yticks(fontsize=14, family='serif')
plt.xlabel('x (km)', fontdict = AxesFont)
plt.ylabel('y (km)', fontdict = AxesFont)
plt.title( title,    fontdict = TitleFont)
plt.axis('equal')
plt.savefig("figure3_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plots 4 & 5: femtosat position error ###########

title    = 'Femtosatellite Position Estimate Error'
y_label  = 'Error (km)'
# Residual  = estimate - truth
x_resids = x_est - x_truth
y_resids = y_est - y_truth
resids   = (x_resids**2 + y_resids**2)**0.5

plt.plot(states[:,0]/60.0, resids, 'k.')
plt.grid('on')
plt.xticks(fontsize=14, family = 'serif')
plt.yticks(fontsize=14, family = 'serif')
plt.xlabel(x_label,     fontdict = AxesFont)
plt.ylabel(y_label,     fontdict = AxesFont)
plt.title(title,        fontdict = TitleFont)
plt.savefig("figure4_case" + str(test_case) + ".png")
plt.show()
fig.clear()

# Now we also want to see what the error would have been had we just
# propagated the initial guess, without filtering any measurements.
# So start by re-plotting the same KF estimate residuals.
plt.plot(states[:,0]/60.0, resids, 'k-', label='KF Residuals')
plt.grid('on')
plt.xticks(fontsize=14,  family = 'serif')
plt.yticks(fontsize=14,  family = 'serif')
plt.xlabel(x_label,      fontdict = AxesFont)
plt.ylabel('Error (km)', fontdict = AxesFont)
plt.title(title,         fontdict = TitleFont)
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
initial_error  = ((x0_guess[0]-x0_fem_truth[0])**2 +
                  (x0_guess[1]-x0_fem_truth[1])**2)**0.5
print "True initial position : %.3f" % x0_fem_truth[0], " %.3f" % \
      x0_fem_truth[1]
print "Initial position guess: %.3f" % x0_guess[0], " %.3f" % \
      x0_guess[1]
print "Initial guess error   : %.3f" % initial_error
plt.plot(states[:,0]/60.0, resids_2body, 'r-', 
         label='Guess Propagation Residuals')
plt.legend(prop=LegFont, shadow='True', loc='best')
plt.savefig("figure5_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 6: femtosat velocity error ###########

title       = 'Femtosatellite Velocity Estimate Error'
y_label     = 'Error (m/s)'
# Residual  = estimate - truth
xdot_resids = xd_est - xd_truth
ydot_resids = yd_est - yd_truth
vel_resids  = (xdot_resids**2 + ydot_resids**2)**0.5

plt.plot(states[:,0]/60.0, vel_resids*1000.0,  'k.')
plt.grid('on')
plt.xticks(fontsize=14, family = 'serif')
plt.yticks(fontsize=14, family = 'serif')
plt.xlabel(x_label,     fontdict = AxesFont)
plt.ylabel(y_label,     fontdict = AxesFont)
plt.title(title,        fontdict = TitleFont)
plt.savefig("figure6_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 7: position error & 3sig bounds ###########

# Error  = estimate - truth
time_array = states[:,0]/60.0
fig = plt.figure()
ax1 = fig.add_subplot(211)

# 3*sigma bounds
Three_sig_11_f = 3*(states[:,13])**0.5;
Three_sig_22_f = 3*(states[:,14])**0.5;
Three_sig_11_b = 3*(states[:,17])**0.5;
Three_sig_22_b = 3*(states[:,18])**0.5;

y_label = 'x error (km)'
ax1.plot(time_array, x_resids,       'k-' )
ax1.plot(time_array, Three_sig_11_f, 'g--')
ax1.plot(time_array,-Three_sig_11_f, 'g--')
ax1.plot(time_array, Three_sig_11_b, 'b--')
ax1.plot(time_array,-Three_sig_11_b, 'b--')
plt.xticks(fontsize=14, family = 'serif')
plt.yticks(fontsize=14, family = 'serif')
plt.ylabel(y_label,     fontdict = AxesFont)
plt.title(u'Position Error and 3\u03C3 Bounds: Forward-Backward\ngreen = forward, blue = backwards', fontdict = TitleFont)

y_label = 'y error (km)'
ax2 = fig.add_subplot(212)
ax2.plot(time_array, y_resids,       'k-' )
ax2.plot(time_array, Three_sig_22_f, 'g--')
ax2.plot(time_array,-Three_sig_22_f, 'g--')
ax2.plot(time_array, Three_sig_22_b, 'b--')
ax2.plot(time_array,-Three_sig_22_b, 'b--')
plt.xticks(fontsize=14, family = 'serif')
plt.yticks(fontsize=14, family = 'serif')
plt.xlabel(x_label,     fontdict = AxesFont)
plt.ylabel(y_label,     fontdict = AxesFont)
plt.savefig("figure7_case" + str(test_case) + ".png")
plt.show()
fig.clear()

#################### Plot 8: velocity error & 3sig bounds ###########

fig = plt.figure()
ax1 = fig.add_subplot(211)

# 3*sigma bounds
Three_sig_33_f = 3*(states[:,15])**0.5;
Three_sig_44_f = 3*(states[:,16])**0.5;
Three_sig_33_b = 3*(states[:,19])**0.5;
Three_sig_44_b = 3*(states[:,20])**0.5;

y_label = 'x_dot error (km/s)'
ax1.plot(time_array, xdot_resids,    'k-' )
ax1.plot(time_array, Three_sig_33_f, 'g--')
ax1.plot(time_array,-Three_sig_33_f, 'g--')
ax1.plot(time_array, Three_sig_33_b, 'b--')
ax1.plot(time_array,-Three_sig_33_b, 'b--')
plt.xticks(fontsize=14, family = 'serif')
plt.yticks(fontsize=14, family = 'serif')
plt.ylabel(y_label,     fontdict = AxesFont)
plt.title(u'Velocity Error and 3\u03C3 Bounds: Forward-Backward\ngreen = forward, blue = backwards', fontdict = TitleFont)

y_label = 'y_dot error (km/s)'
ax2 = fig.add_subplot(212)
ax2.plot(time_array, ydot_resids,    'k-' )
ax2.plot(time_array, Three_sig_44_f, 'g--')
ax2.plot(time_array,-Three_sig_44_f, 'g--')
ax2.plot(time_array, Three_sig_44_b, 'b--')
ax2.plot(time_array,-Three_sig_44_b, 'b--')
plt.xticks(fontsize=14, family = 'serif')
plt.yticks(fontsize=14, family = 'serif')
plt.xlabel(x_label,     fontdict = AxesFont)
plt.ylabel(y_label,     fontdict = AxesFont)
plt.savefig("figure8_case" + str(test_case) + ".png")
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
#  14) P_11 error covariance, forward     13
#  15) P_22 error covariance, forward     14
#  16) P_33 error covariance, forward     15
#  17) P_44 error covariance, forward     16
#  18) P_11 error covariance, back pass   17
#  19) P_22 error covariance, back pass   18
#  20) P_33 error covariance, back pass   19
#  21) P_44 error covariance, back pass   20
#
# footnote [2]
#   To plot a blue circle instead of Earth image, use:
#     Earth_outline = plt.Circle((0,0), r_Titan, color='blue')
#     fig.gca().add_artist(Earth_outline)