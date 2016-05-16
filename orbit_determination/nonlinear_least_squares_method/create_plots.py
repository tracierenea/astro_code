#!/usr/bin/env python
# 
# Tracie Perez
#
# Plot the synthetic measurements and estimate results.

from numpy import *
import matplotlib.pyplot as plt

# Read in data from file
data   = loadtxt('measurement_data.txt') # time, y_true, y_meas
states = loadtxt('sat_states.txt') # size m x 13

# Create plot of measurement data
fig = plt.gcf()
plt.plot(data[:,0]/60.0, data[:,2], 'o', label='Measurements', 
  alpha=0.6, markerfacecolor='blue', markeredgecolor='blue')
plt.plot(data[:,0]/60.0, data[:,1], 'k-', label='Truth')
plt.xlabel('Time (minutes)', fontsize=14)
plt.ylabel('Doppler Shift (Hz)', fontsize=14)
plt.legend(shadow=True, fontsize=14)
# plt.savefig("case3_MeasurementData.png")
plt.show()
fig.clear()


# Create plots of satellite positions (truth and estimates)
r_Earth = 6378.0  # km, radius of Earth
Earth_outline = plt.Circle((0,0), r_Earth, color='blue')
fig = plt.gcf()
fig.gca().add_artist(Earth_outline)

plt.plot(states[:,1], states[:,2], 'ko', label='Mothersat')

plt.plot(states[:,5], states[:,6], 's', label='Femtosat (truth)',
         markerfacecolor='green', markeredgecolor='green')

plt.plot(states[:,9], states[:,10], 's',
         label='Femtosat (NLS estimate)', markerfacecolor='red',
         markeredgecolor='red')

plt.axis('equal')
plt.grid('on')
plt.legend(shadow=True, fontsize=14, loc='best')
plt.xlabel('x (km)', fontsize=14)
plt.xlabel('y (km)', fontsize=14)
# plt.savefig("case3_StateData.png")
plt.show()