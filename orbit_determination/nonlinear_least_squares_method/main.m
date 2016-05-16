% Tracie Perez
%
% Orbit determination from Doppler shift measurements.
%
% The purpose of this simulation is to estimate the initial state of
% the femtosatellite's orbit (r, r_dot) using Doppler shift
% measurements made on the mothersat.
%
% These scripts are configured for octave, but can be used with
% MATLAB if adjustments are made to the ode45 calls.
%
% My preferred way to run this on the terminal command line is:
% > clear; octave main.m; python create_plots.py 
%

clc; close all; clear all;

% Octave package odepkg must be installed and loaded here
pkg load odepkg

% Ignore time stamps for all function files, this improves speed
ignore_function_time_stamp = "all";

% Constants chosen for this simulation
mu          = 398600;                % km^3/s^2, grav param. of Earth
freq        = 435e6;                 % Hz, frequency of tx'ed signal
rad_mom     = 26600;                 % km, the GPS S/V orbit
deploy_v    = 2/1000;                % km/sec (deployment velocity)
t0          = 0;                     % initial time
dt          = 1;                     % seconds, time interval
m           = 1000;                   % how many measurements to make
guess_error = [100; 100; 2; 2];    % add to true IC to given guess
noise_std   = 10;                     % standard deviation of noise

% Initial state for mothersat
y_dot0_mom  = sqrt(mu/rad_mom);      % km/sec (circular orbit)
X0_mom      = [rad_mom; 0; 0; y_dot0_mom];

% A few different initial conditions for the femtosat are considered.
% Uncomment one of these cases to be analyzed.
%------------------------------------------------------------------%

% Test 1: right after deployment from mothersat
% test       = 1;
% rad_fem    = rad_mom;
% y_dot0_fem = sqrt(mu/rad_fem);  % km/sec (circular orbit)
% X0_fem     = [rad_fem; 0; -deploy_v; y_dot0_fem];

% Test 2: femsat already decayed in altitude, but on x-axis
% test       = 2;
% rad_fem    = 400 + 6371;        % km, the ISS orbit
% y_dot0_fem = sqrt(mu/rad_fem);  % km/sec (circular orbit)
% X0_fem     = [rad_fem; 0; -deploy_v; y_dot0_fem];

% Test 3: mom at 30 degrees and fem at 60 degrees off of x-axis (not
% including the deployment velocity component)
test       = 3;
rad_fem    = 400 + 6378;        % km, the ISS orbit
y_dot0_fem = sqrt(mu/rad_fem);  % km/sec (circular orbit)
vel_fem    = y_dot0_fem;
vel_mom    = y_dot0_mom;
X0_fem     = [ rad_fem*cosd(60);
               rad_fem*sind(60);
              -vel_fem*cosd(30);
               vel_fem*sind(30)];
X0_mom     = [ rad_mom*cosd(30);
               rad_mom*sind(30);
              -vel_mom*cosd(60);
               vel_mom*sind(60)];

%------------------------------------------------------------------%

% Initial state guess: [r_x r_y r_x_dot r_y_dot]'
x0_true         = X0_fem;
x0_guess        = X0_fem + guess_error;

% Generate truth data
tf   = t0 + (m-1)*dt;                 % seconds, final time
time = [t0:dt:tf]';                   % time vector for analysis
[~, fem_states] = ode45(@two_body_EOM, time, X0_fem, mu);
[~, mom_states] = ode45(@two_body_EOM, time, X0_mom, mu);

% For the test case #1 scenario: at the deployment epoch, the range
% between the sats is 0, and we'll have a divide by 0 error. So,
% remove the NaN measurement data point of the first epoch
if test == 1
    mom_states = mom_states(2:end,:);
    fem_states = fem_states(2:end,:);
    time       = time(2:end,:);
    m          = m-1;
end

% Generate measurement data
y_true     = get_measurements(fem_states, mom_states, freq);
y_meas     = y_true + noise_std*randn(size(y_true));

% Initialize matrices and vectors
n          = length(x0_guess);        % Number of states
Phi0_vec   = reshape(eye(n), n^2, 1);
max_count  = 10;                      % Max # of iterations to allow
R_big      = noise_std^2*eye(m);      % Obtain R (covariance) matrix 
invR       = inv(R_big);              % Compute inverse of R matrix

% Print out some info for confirmation
fprintf('Read in %i measurements.\n',m);
fprintf('Std of noise specified to be: %.1f Hz\n',noise_std);
fprintf('Initial state of mothersat: [%7.1f %7.1f %7.3f %7.3f]\n',...
        X0_mom);
fprintf('Initial state of femsat:    [%7.1f %7.1f %7.3f %7.3f]\n',...
        X0_fem);
fprintf('Initial state guess:        [%7.1f %7.1f %7.3f %7.3f]\n',...
        x0_guess);
fprintf('Initial state guess error:  [%7.1f %7.1f %7.3f %7.3f]\n',...
        guess_error);
disp('')
data_out = zeros(m,13);
data_out(:,1) = time';
data_out(:,2:5) = mom_states;
data_out(:,6:9) = fem_states;

% Call the main NLS function to solve for initial state estimate
[x_estimate, fem_state_estimates] = ...
  nonlinear_least_squares_orbit_det(time, x0_guess, y_meas,Phi0_vec, invR, max_count, mom_states, freq, mu);

print_results(x0_true, x0_guess, x_estimate);
data_out(:,10:13) = fem_state_estimates;

% Save data for later plotting (I prefer using python's matplotlib
% to octave's gnuplot)
file_id = fopen('measurement_data.txt', 'w');
for ii = 1:m
  fprintf(file_id,'%.6e %.6e %.6e\n',time(ii),y_true(ii),y_meas(ii));
end
file_id = fopen('sat_states.txt','w');
for ii = 1:m
  fprintf(file_id,'%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n', data_out(ii,:));
end
status = fclose("all");