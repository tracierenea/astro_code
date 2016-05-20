              Orbit determination from Doppler measurements

Assumptions:
- Planetary body and satellites are assumed to be point-masses (no
  interruption in the line-of-sight, i.e. the signal, between them)
- Model is just 2-body equation of motion; no other forces considered
- Goal is to estimate the initial condition (position and velocity)
  of the femtosatellite. No other parameters being estimated at this
  time.
- Experimenting with simulation parameters such as: sample rate,
  initial guess, standard deviation of the noise, and convergence 
  criteria.
  
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
- Because I prefer the flexibility of pyplot over gnuplot, plots are
  created in the python file create_plots.py

My preferred way to run this on the terminal command line is:
  > clear; octave main.m; python create_plots.py 
