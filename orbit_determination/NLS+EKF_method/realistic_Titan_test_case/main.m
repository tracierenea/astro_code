function [] = main(test_case)
%MAIN Orbit determination from Doppler shift measurements.
%
%  Author: Tracie Perez
%
%  The purpose of this simulation is to estimate the initial state
%  of the femtosatellite's orbit (r, r_dot) using Doppler shift
%  measurements made on the mothersat.
%
%  One argument required:
%    test_case : 1, 2, or 3
%
%  These scripts are configured for octave, but can be used with
%  MATLAB if adjustments are made to the ode45 calls.
%
%

% Constants chosen for this simulation

G          = 6.6742e-11 / 1000**3;% km^3/kg*s^2, Univ. grav. constant
m_Titan    = 1.3452e23;           % kg, mass of Titan
mu         = G*m_Titan;           % km^3/s^2, Titan's grav. parameter
rad_Titan  = 2575.5;              % km, radius of Titan
% freq       = 435e6;               % Hz, freq of transmitted signal
freq = 5.65e9; % 5.650 GHz
rad_mom    = 1500+rad_Titan;      % km
deploy_v   = 2/1000;              % km/sec (deployment velocity)
t0         = 0;                   % initial time
dt         = 1;                   % seconds, time interval
m          = 10000;               % try to make this many meas's
max_count  = 3;                   % Max # of iterations to allow
m_NLS      = 1000;                % meas's for NLS

% Standard deviation of noise
if     test_case == 1
  noise_std  = 10;
elseif test_case == 2
  noise_std  = 500;
elseif test_case == 3
  noise_std  = 10;
end
%-------------------------------------------------------------------%

% Octave package odepkg must be installed and loaded here
pkg load odepkg

% Ignore time stamps for all function files, this improves speed
ignore_function_time_stamp = "all";

% Initial state for mothersat
y_dot0_mom = sqrt(mu/rad_mom);      % km/sec (circular orbit)
X0_mom     = [rad_mom; 0; 0; y_dot0_mom];

% A few different initial conditions for the femtosat are considered.
if test_case == 1
  % Test 1: right after deployment from mothersat (thus, error on 
  % the guess isn't as high as other cases...)
  rad_fem     = rad_mom;
  y_dot0_fem  = sqrt(mu/rad_fem);    % km/sec (circular orbit)
  X0_fem      = [rad_fem; 0; -deploy_v; y_dot0_fem];
  guess_error = [25*randn; 25*randn; .1*randn; .1*randn];
elseif test_case == 2
  % Test 2: femsat already decayed in altitude, position on x-axis
  rad_fem     = 400 + rad_Titan;     % km
  y_dot0_fem  = sqrt(mu/rad_fem);    % km/sec (circular orbit)
  X0_fem      = [rad_fem; 0; 0; y_dot0_fem];
  guess_error = [25*randn; 25*randn; 0.1*randn; 0.1*randn];
elseif test_case == 3
  % Test 3: mom at 30 degrees and fem at 60 degrees off of x-axis
  % (not including the deployment velocity component)
  rad_fem     = 400 + rad_Titan;     % km
  y_dot0_fem  = sqrt(mu/rad_fem);    % km/sec (circular orbit)
  vel_fem     = y_dot0_fem;
  vel_mom     = y_dot0_mom;
  X0_fem      = [ rad_fem*cosd(60);
                  rad_fem*sind(60);
                 -vel_fem*cosd(30);
                  vel_fem*sind(30)];
  X0_mom      = [ rad_mom*cosd(30);
                  rad_mom*sind(30);
                 -vel_mom*cosd(60);
                  vel_mom*sind(60)];
  guess_error = [25*randn; 25*randn; 0.1*randn; 0.1*randn];

end

% Initial state guess: [r_x r_y r_x_dot r_y_dot]'
x0_true  = X0_fem;
x0_guess = X0_fem + guess_error;

% Generate truth data
m        = m+1;                          % b/c we'll delete first one
tf       = t0 + (m-1)*dt;                % seconds, final time
time_vec = [t0:dt:tf]';                  % time vector for analysis
[~, fem_states] = ode45(@two_body_EOM, time_vec, X0_fem, mu);
[~, mom_states] = ode45(@two_body_EOM, time_vec, X0_mom, mu);

% For the test case #1 scenario: at the deployment epoch, the range
% between the sats is 0, and we'll have a divide by 0 error. So,
% remove the NaN measurement data point of the first epoch
if test_case == 1
    mom_states = mom_states(2:end,:);
    fem_states = fem_states(2:end,:);
    time_vec   = time_vec(2:end,:);
    m          = m-1;
end

% Generate measurement data
meas_data = get_measurements(time_vec, fem_states, mom_states, ...
                             freq, rad_Titan, 1);
time_meas = meas_data(:,1);
y_true    = meas_data(:,2);
y_meas    = y_true + noise_std*randn(size(y_true));

% Initialize matrices and vectors
n        = length(x0_guess);        % Number of states
Phi0_vec = reshape(eye(n), n^2, 1);
R        = noise_std^2;             % Measurement covariance
R_big    = noise_std^2*eye(m_NLS);  % Obtain R (covariance) matrix 
invR     = inv(R_big);              % Compute inverse of R matrix

% Print out some info for confirmation
disp('')
fprintf('\tIterated extended Kalman filter for orbit determination');
fprintf(', test case: %i\n\n', test_case)
fprintf('%i time points evaluated, %i measurements created', m, ...
        size(y_meas,1));
fprintf(' (%.1f %%)\n', (size(y_meas,1)/m)*100);
fprintf('First %i of the total %i measurements to be used for NLS\n',
        m_NLS, size(y_meas,1));
fprintf('dt = %.3f seconds.\n',dt);
fprintf('Std of noise specified to be: %.1f Hz\n',noise_std);
fprintf('Number of iterations for NLS and EKF: %i\n\n', max_count);
fprintf('Initial state of mothersat: [%7.1f %7.1f %7.3f %7.3f]', ...
        X0_mom);
fprintf('  (truth)\n');
fprintf('Initial state of femsat:    [%7.1f %7.1f %7.3f %7.3f]', ...
        X0_fem);
fprintf('  (truth)\n');
fprintf('Initial state guess:        [%7.1f %7.1f %7.3f %7.3f]\n',...
        x0_guess);
fprintf('Initial state guess error:  [%7.1f %7.1f %7.3f %7.3f]\n',...
        guess_error);

% Step 1: Call the NLS routine for the first m_NLS measurements to
% generate an initial state guess for the EKF.

% First - check to see if there are actually measurements for the
% first m_NLS measurements
continue_flag = 1;
if m_NLS > size(y_meas,1)
  fprintf('Error: trying to use NLS on first %i measurements, but there are only %i measurements total.\n', m_NLS, size(y_meas,1));
  continue_flag = 0;
end
if ~continue_flag
  fprintf('Not able to run NLS on first %i measurements. Exiting.\n', m_NLS);
end
if continue_flag
  for index = 2:m_NLS
    if (time_meas(index) - (time_meas(index-1)+dt)) > 1e-12
      fprintf('Error: There is a gap in the first');
      fprintf(' %i measurements.\n', m_NLS);
      continue_flag = 0;
    end
  end
  if ~continue_flag
    fprintf('Trying to use first %i measurements for NLS estimate, but they are not continuous. Try a smaller number. Exiting. \n', m_NLS);
  end
end

% Second - if so, proceed with NLS run for m_NLS measurements
if continue_flag
  [x_est_NLS, fem_state_est_NLS] =                               ...
  nonlinear_least_squares_orbit_det(time_vec(1:m_NLS), x0_guess, ...
  y_meas(1:m_NLS), Phi0_vec, invR, max_count,                    ...
  mom_states(1:m_NLS,:), freq, mu, rad_Titan);
  
  fprintf('\nResult from NLS, to use as guess for EKF:')
  print_results(x0_true, x0_guess, x_est_NLS);

  % Step 2: Use the resulting initial state guess from the NLS method
  % to initialize the EKF algorithm (uses all measurements).
  [x_estimate, fem_state_est_EKF, P_for, P_back] =               ...
  iterated_EKF_orbit_det(time_vec, x_est_NLS, meas_data, R,      ...
  max_count, mom_states, freq, mu, rad_Titan);

  % Print final results
  fprintf('Results from EKF:')
  print_results(x0_true, x_est_NLS, x_estimate);

  % Save data for later plotting. (I prefer using python's matplotlib
  % to octave's gnuplot.)
  file_id = fopen('measurement_data.txt', 'w');
  fprintf(file_id,'%.6e %.6e %.6e %.6e\n', x0_guess);
  for ii = 1:size(meas_data,1)
    fprintf(file_id,'%.6e %.6e %.6e\n',time_meas(ii), y_true(ii), ...
                                       y_meas(ii));
  end

  file_id = fopen('sat_states.txt','w');
  % first line is true initial state at time 0
  fprintf(file_id,'%.6e %.6e %.6e %.6e\n', x0_true);
  data_out          = zeros(m,17);
  data_out(:,1)     = time_vec';            % time of states
  data_out(:,2:5)   = mom_states;           % truth
  data_out(:,6:9)   = fem_states;           % truth
  data_out(:,10:13) = fem_state_est_EKF;    % estimate of fem state
  data_out(:,14)    = P_for(:, 1);          % P_xx covariance forward
  data_out(:,15)    = P_for(:, 4);          % P_yy covariance forward
  data_out(:,16)    = P_back(:,1);          % P_xx covariance backw.
  data_out(:,17)    = P_back(:,4);          % P_yy covariance backw.
  for ii = 1:m
    fprintf(file_id,'%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n', data_out(ii,:));
  end
  status = fclose("all");
else
  disp('Process halted.')
end
