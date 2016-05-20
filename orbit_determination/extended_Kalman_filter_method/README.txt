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
- Measurements are created and extended Kalman filter method is
  applied in the file main.m. This is currently configured for
  octave, but can be ran in Matlab if the ode45 function calls are
  updated accordingly. i.e. change this:
    [~, fem_states] = ode45(@TwoBodyEOM, time, X0_fem, mu);  
  to this:
    [~, fem_states] = ode45(@TwoBodyEOM, time, X0_fem, [], mu);
- Because I prefer the flexibility of pyplot over gnuplot, plots are
  created in the python file create_plots.py

My preferred way to run this on the terminal command line is:
  > clear; octave main.m; python create_plots.py 
