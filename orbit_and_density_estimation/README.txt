Orbit Determination and Atmospheric Density Estimation from Doppler Measurements

Monte Carlo runs aren't setup yet. To kick off a single run: 

  > bash start_sim.sh $opt

where:
- start_sim does one run for the test case conditions hard-coded in main.m
- $opt is the test case to evaluate: 1, 2, or 3

--------------------------------------------------------------------------------

Assumptions:
- Satellites are in orbit around Titan.
- Satellites are considered to be point masses. The model, for now, is just the
  two-body equation of motion; no other forces considered.
- Titan blocking the line-of-sight between mothersat and femtosat prevents
  Doppler measurement (i.e. if no signal can be received by mothersat, then no
  measurement is available to Kalman filter)
- The goal is to estimate the initial condition state (position and velocity)
  of the femtosatellite, as well as the atmospheric density over the duration
  of the simulation.
- Can experiment with simulation parameters such as: sample rate, initial guess,
  standard deviation of the noise, and convergence criteria.

  
Files:
- Measurements are created and the estimation method is completed in
  the script main.m. The estimation method is two-part:
  1) read in a small subset of the measurements and process with NLS
     algorithm to improve the initial state guess
  2) use this improved initial state guess to initialize the EKF,
     which processes all the measurements and estimates the states
     and IC of the femtosatellite

  The script main.m is is currently configured for octave, but can
  be ran in Matlab if the ode45 function calls are updated
  accordingly. i.e. change this:
    [~, fem_states] = ode45(@TwoBodyEOM, time, X0_fem, mu);  
  to this:
    [~, fem_states] = ode45(@TwoBodyEOM, time, X0_fem, [], mu);
- To install octave's ode package (a package for solving ordinary differential
  equations), run this on the octave command line (NOT IN UBUNTU'S TERMINAL):
  pkg install -forge odepkg
- Because I prefer the flexibility of pyplot over gnuplot, plots are
  created in the python file create_plots.py (or create_plots_MC.py)
  
  Preliminary results from these scripts were presented at the 2016
  AIAA/USU Small Satellite Conference in Logan, Utah.


