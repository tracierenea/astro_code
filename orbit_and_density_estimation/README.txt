Orbit Determination and Atmospheric Density Estimation from Doppler Measurements

Estimate a model of the atmospheric density of Titan. Different test cases use
different measurement types.

Directories:

methodology_validation              : uses the full position and velocity state
                                      of the femtosat
                                      
EKF_with_Doppler_measurement        : uses only Doppler shift measurements

EKF_with_Doppler_range_measurements : uses Doppler shift and range (between
                                      mothersat and femtosat) measurements
--------------------------------------------------------------------------------

Assumptions:
- Satellites are in orbit around Titan.
- The dynamic model for the mothersat is as a point mass undergoing the two-body
  equation of motion; no other forces considered.
- The femtosat is model as a tiny spherical satellite experiencing the
  gravitational force of Titan, as well as induced drag from Titan's
  atmosphere.
- It is assumed that measurements are available every second, despite whether or
  not Titan is blocking the line-of-sight between the two satellites.
- The goal is to estimate the state of the femtosatellite (position and 
  velocity) over the duration of the simulation, as well as two constants
  describing the atmospheric density model.
- Can experiment with simulation parameters such as: sample rate, initial guess,
  standard deviation of the noise, and convergence criteria.

These scripts are configured for MATLAB. The Kalman filter algorithm in general
seems to execute 7 times faster in MATLAB than in Octave, so Octave processing
was abandoned.
