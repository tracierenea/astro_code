function [y] = get_measurements(fem_states, mom_states, freq)
%GETMEASUREMENTS Generate the synthetic Doppler shift measurements
%   Given the states of the mothersat and femtosat, generate
%   synthetic Doppler shift measurements.

m     = size(fem_states, 1);          % number of measurements
n     = size(fem_states, 2);          % number of states
c     = 3e8;                          % speed of light, m/s
lamda = (c/freq)/1000;                % wavelength, km

for ii = 1:m
  r            =  [fem_states(ii,1); fem_states(ii,2)];
  R            =  [mom_states(ii,1); mom_states(ii,2)];
  r_dot        =  [fem_states(ii,3); fem_states(ii,4)];
  R_dot        =  [mom_states(ii,3); mom_states(ii,4)];
  rho_dot_mag  =  get_rho_dot_mag(r,r_dot,R,R_dot);
  y(ii,1)      = -rho_dot_mag/lamda; % Doppler shift
end