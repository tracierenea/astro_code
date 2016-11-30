function [y] = get_measurements(time, fem_states, mom_states, freq)
%GETMEASUREMENTS Generate the synthetic Doppler shift measurements
%   Given the states of the mothersat and femtosat, generate
%   synthetic Doppler shift measurements.
% 
%   Inputs:
%     time          : array of time instances for fem and mom states
%     fem_states    : matrix of the femtosatellite states
%     mom_states    : matrix of the mothersatellite's states
%     freq          : transmission frequency
%
%   Output:
%     y             : 3-column matrix, [time, Doppler shift, range]
%

m             = size(fem_states, 1);         % number of measurements
n             = size(fem_states, 2);         % number of states
c             = 3e8;                         % speed of light, m/s
lamda         = (c/freq)/1000;               % wavelength, km
y             = zeros(0,3);                  % create empty array
meas_index    = 1;

for ii = 1:m
    x_fem       =  fem_states(ii,1);
    y_fem       =  fem_states(ii,2);
    z_fem       =  fem_states(ii,3);
    x_mom       =  mom_states(ii,1);
    y_mom       =  mom_states(ii,2);
    z_mom       =  mom_states(ii,3);
    r           =  fem_states(ii,1:n/2)';
    R           =  mom_states(ii,1:n/2)';
    r_dot       =  fem_states(ii,n/2+1:n)';
    R_dot       =  mom_states(ii,n/2+1:n)';
    rho         =  r - R;                    % Range between mother and fem
    rho_mag     =  norm(rho);
    rho_dot_mag =  (1/rho_mag)*(r'*r_dot - r'*R_dot - R'*r_dot + R'*R_dot);
    Doppler     = -rho_dot_mag/lamda;        % Doppler shift
    
    y(meas_index, :) = [time(ii) Doppler rho_mag]; 
    meas_index       = meas_index + 1;

end % for ii = 1:m
