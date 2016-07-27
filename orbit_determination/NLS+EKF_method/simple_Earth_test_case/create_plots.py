#!/usr/bin/env python
# 
# Tracie Perez
#
# Plot the synthetic measurements and state estimate results.
#

from numpy import *
import matplotlib.pyplot as plt
import os

#### Read in data from file
meas_data   = loadtxt('measurement_data.txt')  # time, y_true, y_meas
states      = loadtxt('sat_states.txt')        # size m x 13. see [1]

#### General parameters/settings
r_Earth = 6378.0  # km, radius of Earth
 # The full path to the image file must be specified for the function
 # get_sample_data, otherwise it assumes the path is relative to the
 # mpl-data/sample_data directory
Earth_image_file = os.getcwd() + \
                   "/The_Earth_seen_from_Apollo_17_cropped.jpg"
# Read in image and plot Earth
image = plt.imread(Earth_image_file)

#### Create plot of measurement data
fig = plt.gcf()
plt.plot(meas_data[:,0]/60.0, meas_data[:,2], 'o', 
         label='Measurements', alpha=0.6, markerfacecolor='blue',
         markeredgecolor='blue')
plt.plot(meas_data[:,0]/60.0, meas_data[:,1], 'k-', label='Truth')
plt.xlabel('Time (minutes)', fontsize=14)
plt.ylabel('Doppler Shift (Hz)', fontsize=14)
plt.title('Measurements Created for Simulation', fontsize=16)
plt.legend(shadow=True, fontsize=14, loc='best')
# plt.savefig("filename.png")
plt.show()
fig.clear()

#### Create plots of satellite positions (truth and estimates)
fig = plt.gcf()
# Add image of Earth
plt.imshow(image, extent=[-r_Earth,r_Earth,-r_Earth,r_Earth])
# Plot all points
plt.plot(states[:,1], states[:,2], 'ko', label='Mothersat')
plt.plot(states[:,5], states[:,6], 's', label='Femtosat (truth)',
         markerfacecolor='green', markeredgecolor='green')
plt.plot(states[:,9], states[:,10], 's',
         label='Femtosat (NLS+EKF estimate)', markerfacecolor='red',
         markeredgecolor='red')
plt.legend(shadow=True, fontsize=14, loc='upper left')
# Start and end points
end = states.shape[0] - 1
plt.text(states[0,1]-1500,   states[0,2]-1500,  "Start")
plt.text(states[end,1]-1500, states[end,2]+500, "End")
plt.text(states[0,9]+500,   states[0,10],  "Start")
plt.text(states[end,9]+500, states[end,10], "End")
plt.axis('equal')
plt.xlabel('x (km)', fontsize=14)
plt.xlabel('y (km)', fontsize=14)
plt.title('Estimated and True Femtosatellite Position', fontsize=16)
plt.show()
fig.clear()

#### Create plot showing when measurements are available
plt.imshow(image, extent=[-r_Earth,r_Earth,-r_Earth,r_Earth])
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
plt.text(states[0,1]-1500,   states[0,2]-1500,  "Start")
plt.text(states[end,1]-1500, states[end,2]+500, "End")
plt.text(states[0,9]+500,   states[0,10],  "Start")
plt.text(states[end,9]+500, states[end,10], "End")   
plt.axis('equal')
plt.xlabel('x (km)', fontsize=14)
plt.xlabel('y (km)', fontsize=14)
plt.title('Measurement Availability\nGreen = possible, Red = not possible due to signal blockage', fontsize=16)
plt.show()
fig.clear()

#### Create plot of femtosat position error (truth - estimate)
x_resids = states[:,5] - states[:,9]
y_resids = states[:,6] - states[:,10]
resids   = (x_resids**2 + y_resids**2)**0.5
plt.plot(states[:,0]/60.0, resids, 'k.')
plt.grid('on')
plt.xlabel('Time (minutes)', fontsize=14)
plt.ylabel('Error (km)', fontsize=14)
plt.title('Femtosatellite Position Estimate Error', fontsize=16)
plt.show()



# footnote [1]
#   The columns of the array states will be, from left to right:
#   1) time of the state
#   2) r_x     of mothersat (truth)
#   3) r_y     of mothersat (truth)
#   4) r_x_dot of mothersat (truth)
#   5) r_y_dot of mothersat (truth)
#   6) r_x     of femsat    (truth)
#   7) r_y     of femsat    (truth)
#   8) r_x_dot of femsat    (truth)
#   9) r_y_dot of femsat    (truth)
#  10) r_x     of femsat    (estimate)
#  11) r_y     of femsat    (estimate)
#  12) r_x_dot of femsat    (estimate)
#  13) r_y_dot of femsat    (estimate)
#
# footnote [2]
#   To plot a blue circle instead of Earth image, use:
#     Earth_outline = plt.Circle((0,0), r_Earth, color='blue')
#     fig.gca().add_artist(Earth_outline)