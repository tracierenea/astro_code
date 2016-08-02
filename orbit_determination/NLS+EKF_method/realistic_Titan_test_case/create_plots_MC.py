#!/usr/bin/env python
# 
# Tracie Perez
#
# Evaluate the Monte Carlo results.
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

# Establish the true initial position (this is from main_MC.m)
G         = 6.6742e-11 / 1000**3;
m_Titan   = 1.3452e23;
mu        = G*m_Titan;
rad_Titan = 2575.5;
rad_mom   = 1500+rad_Titan;
if test_case == 1:
  rad_fem     = rad_mom;
  y_dot0_fem  = sqrt(mu/rad_fem);
  X0_true     = array([rad_fem, 0])
elif test_case == 2:
  rad_fem     = 400 + rad_Titan;
  y_dot0_fem  = sqrt(mu/rad_fem);
  X0_true     = array([rad_fem, 0])
elif test_case == 3:
  rad_fem     = 400 + rad_Titan;
  X0_true     = [ rad_fem*cos(radians(60)),
                  rad_fem*sin(radians(60))]
else:
  exit('Invalid test case argument given. Options are: 1, 2, 3')

#### For plotting Titan image:
r_Titan = 2575.5  # km, radius of Titan
Titan_image_file = os.getcwd() + "/Titan_cropped.png"
image = plt.imread(Titan_image_file)

#### Read in data from file. The columns, with 4 states each, are:
####  [counter  initial_x0_guess  x0_result_NLS  x0_result_final_EKF]
file_name = 'MC_results_case' + str(test_case) + '.txt'
data      = loadtxt(file_name)
x0_guess  = data[:,1]  # r_x(0) - initial guess
y0_guess  = data[:,2]  # r_y(0) - initial guess
print "Opened file: ", file_name
print "Read in", data.shape[0], "points"


#### Throw out outlier solutions: if x(0) or y(0) has magnitude
#    500 km or more greater than the initial guess, toss
num_points   = data.shape[0]
deleted_rows = 0;
for counter in range(0,num_points):
  index  = counter - deleted_rows
  this_x0_est = data[index, 9]
  this_y0_est = data[index, 10]
  if abs(this_x0_est - x0_guess[index]) > 500:
    data = delete(data, index, axis=0)
    deleted_rows = deleted_rows + 1
  elif abs(this_y0_est - y0_guess[index]) > 500:
    data = delete(data, s_[index], axis=0)
    deleted_rows = deleted_rows + 1
print "Deleted %i outlier points." % deleted_rows

#### Redefine data to plot now that outliers have been removed
x0_guess  = data[:,1]  # r_x(0)     - initial guess
y0_guess  = data[:,2]  # r_y(0)     - initial guess
xd0_guess = data[:,3]  # r_x_dot(0) - initial guess
yd0_guess = data[:,4]  # r_y_dot(0) - initial guess
x0_est    = data[:,9]  # r_x(0)     - final estimate
y0_est    = data[:,10] # r_y(0)     - final estimate
xd0_est   = data[:,11] # r_x_dot(0) - final estimate
yd0_est   = data[:,12] # r_y_dot(0) - final estimate


#### Plot the dispersion of the initial state guesses
fig = plt.gcf()
plt.imshow(image, extent=[-r_Titan,r_Titan,-r_Titan,r_Titan])
plt.plot(x0_guess,   y0_guess,  'ko', label='Guesses')
plt.plot(X0_true[0], X0_true[1],'ro', label="Truth")
plt.legend(shadow=True, fontsize=16, loc='upper left')
plt.xlabel('x (km)', fontsize=16)
plt.ylabel('y (km)', fontsize=16)
plt.title('Dispersion of Initial State (Position) Guesses',
           fontsize=22)
plt.axis('equal')
plt.savefig("figure1_case" + str(test_case) + "_MC.png")
plt.show()
fig.clear()

#### Plot the dispersion of the initial state estimates
plt.imshow(image, extent=[-r_Titan,r_Titan,-r_Titan,r_Titan])
plt.plot(x0_est,     y0_est,    'ko', label='Estimates')
plt.plot(X0_true[0], X0_true[1],'ro', label="Truth")
plt.legend(shadow=True, fontsize=16, loc='upper left')
plt.xlabel('x (km)', fontsize=16)
plt.ylabel('y (km)', fontsize=16)
plt.title('Dispersion of Initial State (Position) Estimates',
           fontsize=22)
plt.axis('equal')
plt.savefig("figure2_case" + str(test_case) + "_MC.png")
plt.show()
fig.clear()

#### Plot truth, guesses, and estimates all in one combo plot
plt.plot(x0_guess,   y0_guess,  'ro', label='Guesses')
plt.plot(x0_est,     y0_est,    'go', label='Estimates')
plt.plot(X0_true[0], X0_true[1],'ko', label="Truth")
plt.legend(shadow=True, fontsize=16, loc='upper left')
plt.xlabel('x (km)', fontsize=16)
plt.ylabel('y (km)', fontsize=16)
plt.axis('equal')
plt.savefig("figure3_case" + str(test_case) + "_MC.png")
plt.show()
fig.clear()

#### Print out stats
print "x(0)\ttruth: %.3f"  % X0_true[0],
print "   mean: %.3f"      % x0_est.mean(),
print "   diff: %.3f"      % float(X0_true[0]-x0_est.mean()),
print "   std: %.3f"       % x0_est.std()
print "y(0)\ttruth: %.3f"  % X0_true[1],
print "   mean: %.3f"      % y0_est.mean(),
print "   diff: %.3f"      % float(X0_true[1]-y0_est.mean()),
print "   std: %.3f\n"     % y0_est.std()