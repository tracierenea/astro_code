%  Main simulalation script: orbit and atmospheric density estimation
%  from Doppler shift measurements.
%
%  Author: Tracie Perez // UTA MAE PhD Student
%
%  The purpose of this simulation is to estimate the initial state
%  of the femtosatellite's orbit (r, r_dot) using Doppler shift
%  measurements made on the mothersat.
%
%  This script differs from the 2D case by giving an initial velocity
%  and/or position in the 3rd dimension (i.e. in the z-axis). Also, where
%  the 2D case does a check for measurement availability by evaluating 
%  line-of-sight, that is non-trivial in 3D and has not been implemented
%  yet here.
%
%  Note: The .m files in this folder are configured for MATLAB, not octave.
%

close all; clear all; clc;

% Set the test case here : 1, 2, or 3 (only 1 works for now)
test_case = 4;

%%% Start timer
t_start = tic;

%%% Configure the simulation by defining the following constants:
G          = 6.6742e-11 / 1000^3; % km^3/kg*s^2, Univ. grav. constant
m_Titan    = 1.3452e23;           % kg, mass of Titan
mu         = G*m_Titan;           % km^3/s^2, Titan's grav. parameter
rad_Titan  = 2575.5;              % km, radius of Titan
freq       = 5.65e9;              % 5.650 GHz, freq of transmitted signal
rad_mom    = 1500+rad_Titan;      % km
deploy_v   = 2/1000;              % km/sec (deployment velocity)
t0         = 0;                   % initial time
dt         = 1;                   % seconds, time interval
m          = 10000;               % try to make this many meas's
max_count  = 3;                   % Max # of iterations to allow
m_NLS      = 1000;                % measurements to use for NLS

%%% Tune the EKF by experimenting with the following parameters:

% Standard deviation of noise and guess error for each case. Note, the
% state vector is defined as: [r_x r_y r_z r_x_dot r_y_dot r_z_dot]'
noise_std =  500;
error_vec = [5; 5; 5; 0.1; 0.1; 0.1];

% Unlike the general EKF, this algorithm isn't sensitive at all to the
% initial covariance matrix provided to the EKF function. This is because
% it process the measurements forwards and backwards in time "max_count"
% number of times, so the effect of the initial covariance does not linger
% around.... Replaced this with 100* and 1000*eye(6) and found the same
% exact estimate of the IC for a given measurement set.
P0 = diag(error_vec.^2);


% A few different initial conditions for the femtosat are considered.
rad_mom    = 1500+rad_Titan;
y_dot0_mom = sqrt(mu/rad_mom);        % km/sec (circular orbit)
if test_case == 1 % right after deployment from mothersat
    X0_mom     = [rad_mom; 0; 0; 0; y_dot0_mom; 0];
    X0_fem     = [rad_mom; 0; 0; 0; y_dot0_mom; deploy_v];
    m          = m+1;                 % b/c we'll delete first one
elseif test_case == 2 % femtosat decayed in altitude, both on x-axis
    X0_mom     = [rad_mom; 0; 0; 0; y_dot0_mom; 0];
    rad_fem    = 400 + rad_Titan;     % km
    y_dot0_fem = sqrt(mu/rad_fem);    % km/sec (circular orbit)
    X0_fem     = [rad_fem; 0; 0; 0; y_dot0_fem; deploy_v];
elseif test_case == 3 % mom @ 60 degrees, femtosat @ 30 deg off x-axis
    rad_fem    = 400 + rad_Titan;     % km
    y_dot0_fem = sqrt(mu/rad_fem);    % km/sec (circular orbit)
    vel_fem    = y_dot0_fem;
    vel_mom    = y_dot0_mom;
    X0_fem     = [ rad_fem*cosd(60);
                   rad_fem*sind(60);
                   0;
                  -vel_fem*cosd(30);
                   vel_fem*sind(30);
                   deploy_v];
    X0_mom     = [ rad_mom*cosd(30);
                   rad_mom*sind(30);
                   0;
                  -vel_mom*cosd(60);
                   vel_mom*sind(60);
                   0];
elseif test_case == 4 % Similar to case 2 but initial position off-plane
    X0_mom     = [rad_mom; 0; 0; 0; y_dot0_mom; 0];
    rad_fem    = 400 + rad_Titan;     % km
    x_dot0_fem = sqrt(mu/rad_fem);    % km/sec (circular orbit)
    X0_fem     = [0; 0; rad_fem; x_dot0_fem; 0; deploy_v];    
               
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial state guess
guess_error     = error_vec.*randn(6,1);
x0_true         = X0_fem;
x0_guess        = x0_true + guess_error;

% Generate truth data
tf              = t0 + (m-1)*dt;             % seconds, final time
time_vec        = [t0:dt:tf]';               % time vector for analysis
[~, mom_states] = ode45(@two_body_EOM, time_vec, X0_mom, [], mu);
[~, fem_states] = ode45(@two_body_EOM, time_vec, X0_fem, [], mu);

% For the test case #1 scenario: at the deployment epoch, the range
% between the sats is 0, and we'll have a divide by 0 error. So,
% remove the NaN measurement data point of the first epoch
if test_case == 1
  mom_states    = mom_states(2:end,:);
  fem_states    = fem_states(2:end,:);
  time_vec      = time_vec(2:end,:);
  m             = m-1;
end

% Generate measurement data
meas_data       = get_measurements(time_vec, fem_states, mom_states, freq);
time_meas       = meas_data(:,1);
y_true          = meas_data(:,2);
y_meas          = y_true + noise_std*randn(size(y_true));

% Initialize matrices and vectors
n               = length(x0_guess);          % Number of states
Phi0_vec        = reshape(eye(n), n^2, 1);   % Reshape matrix into vector
R               = noise_std^2;               % Measurement covariance
R_big           = noise_std^2*eye(m_NLS);    % Obtain R (covariance) matrix
invR            = inv(R_big);                % Compute inverse of R matrix

% Print out some info for confirmation
print_sim_description(test_case, y_meas, m_NLS, error_vec, dt, noise_std,...
                      max_count, guess_error, X0_mom, X0_fem, x0_guess);

% Step 1: Call the NLS routine for the first m_NLS measurements to
% generate an initial state guess for the EKF. First, check to see if there
% are actually measurements for the first m_NLS measurements
continue_flag = 1;
if m_NLS > size(y_meas,1)
  fprintf('Error: trying to use NLS on first ');
  fprintf('%i measurements, but there are only %i measurements total.\n',...
          m_NLS, size(y_meas,1));
  continue_flag = 0;
end
if ~continue_flag
  fprintf('Not able to run NLS on first %i measurements.\n', m_NLS);
  exit;
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
    fprintf('Trying to use first %i measurements for NLS estimate', m_NLS);
    fprintf(', but they are not continuous. Try a smaller number. \n');
    exit;
  end
end

% Second: if so, proceed with NLS run for m_NLS measurements
[x_est_NLS, fem_state_est_NLS] = nonlinear_least_squares_orbit_det(     ...
time_vec(1:m_NLS), x0_guess, y_meas(1:m_NLS), Phi0_vec, invR, max_count,...
mom_states(1:m_NLS,:), freq, mu);

fprintf('\nResult from NLS:\n')
print_results(x0_true, x0_guess, x_est_NLS);

% Step 2: Check the resulting initial state guess from the NLS method.
% If any state estimate seems "way off" go with the original guess instead.
flag_dont_use_NLS_result = 0;
for index = 1:n
    if abs(x_est_NLS(index)) > abs(50*x0_guess(index))
        flag_dont_use_NLS_result = 1;
    end
end
if flag_dont_use_NLS_result
    guess_for_EKF = x0_guess;
    disp('NLS estimate seems worse. Using original IC guess.');
else
    guess_for_EKF = x_est_NLS;
    disp('NLS result passed sanity check; using it.');
end

% Step 3: Run EKF algorithm (uses all measurements)
[x_estimate, fem_state_est_EKF, P_for, P_back] = iterated_EKF_orbit_det(...
time_vec, guess_for_EKF, meas_data, R, max_count, mom_states, freq, mu, P0);

fprintf('Results from F-B EKF:\n')
print_results(x0_true, guess_for_EKF, x_estimate);



% Step 4: propogate the Kalman filter's estimate of the initial state
% forward in time using the two-body equation (not accounting for drag).
[~, fem_states_IC_prop] = ode45(@two_body_EOM, time_vec, x_estimate,[],mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create plots

% Set default color order. This will persist until Matlab is closed.
% Override with individual parameters in plot commands. Using UTA colors.
set(groot,'defaultAxesColorOrder',   ...
    [0       100/255 177/255;        ... % UTA blue
     245/255 128/255  38/255;        ... % UTA orange
     105/255  41/255  23/255;        ... % UTA secondary maroon
     0        68/255 124/255]);      ... % UTA secondary blue

% Plot 1: measurement data: truth and synthetically created meas's
figure(1);
plot(time_vec/60, y_meas, 'o', time_vec/60, y_true, 'k-', 'LineWidth', 2);
xlabel('Time (minutes)');
ylabel('Doppler Shift (Hz)');
set(gca,'FontSize', 16,'XLim',[0 max(time_vec/60)]);
legend('Measurements', 'Truth', 'Location', 'Best');

% Plot 2: satellite positions
figure(2); 
plot3(fem_states(:,1), fem_states(:,2), fem_states(:,3), ...
      'k--', 'LineWidth', 2); hold on;
plot3(fem_state_est_EKF(:,1),  fem_state_est_EKF(:,2),   ...
      fem_state_est_EKF(:,3),  'LineWidth', 2);
plot3(fem_states_IC_prop(:,1), fem_states_IC_prop(:,2),  ...
      fem_states_IC_prop(:,3), 'LineWidth', 2);
% text(fem_state_est_EKF(1,1), fem_state_est_EKF(1,2),  ...
%      fem_state_est_EKF(1,3)-200, 'Start', 'FontSize', 16);
% text(fem_state_est_EKF(end,1), fem_state_est_EKF(end,2),  ...
%      fem_state_est_EKF(end,3)-200, 'End', 'FontSize', 16); 
xlabel('$r_{x}\:\left(km\right)$','Interpreter','LaTex');
ylabel('$r_{y}\:\left(km\right)$','Interpreter','LaTex');
zlabel('$r_{z}\:\left(km\right)$','Interpreter','LaTex');
legend('Truth', 'EKF Estimate', 'IC Estimate Propagation', ...
       'Location', 'Best');
set(gca,'FontSize',16);

% Plot 3: position state estimate errors. 
for index = 1:m
    residuals_EKF(index,:)     = fem_state_est_EKF(index,:) - ...
                                 fem_states(index,:);
    residuals_IC_prop(index,:) = fem_states_IC_prop(index,:) - ...
                                 fem_states(index,:);
    sig3_EKF(index,:)          = 3.*P_back(index,:).^0.5;                   
end
figure(3);
subplot(311); 
plot(time_vec/60, residuals_EKF(:,1),     'LineWidth', 2); hold on;
plot(time_vec/60, residuals_IC_prop(:,1), 'LineWidth', 2);
% plot(time_vec/60, sig3_EKF(:,1), 'k-',    'LineWidth', 2);
% plot(time_vec/60,-sig3_EKF(:,1), 'k-',    'LineWidth', 2);
ylabel('$r_{x}\:\left(km\right)$','Interpreter','LaTex');
legend('EKF Estimate', 'IC Solution Propagation', 'Location', 'Best');
set(gca,'FontSize',16, 'XTickLabel', [], 'XLim',[0 max(time_vec/60)]);
subplot(312); 
plot(time_vec/60, residuals_EKF(:,2),     'LineWidth', 2); hold on;
plot(time_vec/60, residuals_IC_prop(:,2), 'LineWidth', 2);
ylabel('$r_{y}\:\left(km\right)$','Interpreter','LaTex');
set(gca,'FontSize',16, 'XTickLabel', [], 'XLim',[0 max(time_vec/60)]);
subplot(313); 
plot(time_vec/60, residuals_EKF(:,3),     'LineWidth', 2); hold on;
plot(time_vec/60, residuals_IC_prop(:,3), 'LineWidth', 2);
ylabel('$r_{z}\:\left(km\right)$','Interpreter','LaTex');
xlabel('Time (minutes)');
set(gca,'FontSize',16,'XLim',[0 max(time_vec/60)]);

% Plot 4: residuals in velocity state estimates 
figure(4);
subplot(311); 
plot(time_vec/60, residuals_EKF(:,4)*1000,     'LineWidth', 2); hold on;
plot(time_vec/60, residuals_IC_prop(:,4)*1000, 'LineWidth', 2);
% plot(time_vec/60, sig3_EKF(:,4), 'k-',    'LineWidth', 2);
% plot(time_vec/60,-sig3_EKF(:,4), 'k-',    'LineWidth', 2);
ylabel('$\dot{r}_{x}\:\left(\frac{m}{s}\right)$','Interpreter','LaTex');
legend('EKF Estimate', 'IC Solution Propagation', 'Location', 'Best');
set(gca,'FontSize',16, 'XTickLabel', [], 'XLim',[0 max(time_vec/60)]);
subplot(312); 
plot(time_vec/60, residuals_EKF(:,5)*1000,     'LineWidth', 2); hold on;
plot(time_vec/60, residuals_IC_prop(:,5)*1000, 'LineWidth', 2);
ylabel('$\dot{r}_{y}\:\left(\frac{m}{s}\right)$','Interpreter','LaTex');
set(gca,'FontSize',16, 'XTickLabel', [], 'XLim',[0 max(time_vec/60)]);
subplot(313); 
plot(time_vec/60, residuals_EKF(:,6)*1000,     'LineWidth', 2); hold on;
plot(time_vec/60, residuals_IC_prop(:,6)*1000, 'LineWidth', 2);
ylabel('$\dot{r}_{z}\:\left(\frac{m}{s}\right)$','Interpreter','LaTex');
xlabel('Time (minutes)');
set(gca,'FontSize',16, 'XLim',[0 max(time_vec/60)]);


% Plot 5: combo plot with measuremetns and position error
for index = 1:m
    resids_pos_EKF(index)     = norm(residuals_EKF(index,:));
    resids_pos_IC_prop(index) = norm(residuals_IC_prop(index,:)); 
end
figure(5);
yyaxis left;
plot(time_vec/60, y_meas, 'o');
ax = gca;
xlabel('Time (minutes)');
ylabel('Doppler Shift Measurements (Hz)');
yyaxis right
plot(time_vec/60, resids_pos_EKF, 'k', 'LineWidth', 2);
ax = gca;
ax.YColor = 'k';
ylabel('EKF Position Estimate Error (km)');
set(gca,'FontSize',16, 'XLim',[0 max(time_vec/60)]);

% Plot 6: sstimated and true position states
figure(6);
subplot(311); 
plot(time_vec/60, fem_state_est_EKF(:,1),  'LineWidth', 2); hold on;
plot(time_vec/60, fem_states_IC_prop(:,1), 'LineWidth', 2);
plot(time_vec/60, fem_states(:,1),         'LineWidth', 2);
legend('EKF Estimate', 'IC Solution Propagation', 'Truth', ...
       'Location', 'Best');
ylabel('$r_{x}\:\left(km\right)$','Interpreter','LaTex');
set(gca,'FontSize', 16, 'XTickLabel', [], 'XLim',[0 max(time_vec/60)]);
subplot(312);
plot(time_vec/60, fem_state_est_EKF(:,2),  'LineWidth', 2); hold on;
plot(time_vec/60, fem_states_IC_prop(:,2), 'LineWidth', 2);
plot(time_vec/60, fem_states(:,2),         'LineWidth', 2);
ylabel('$r_{y}\:\left(km\right)$','Interpreter','LaTex');
set(gca,'FontSize', 16, 'XTickLabel', [], 'XLim',[0 max(time_vec/60)]);
subplot(313);
plot(time_vec/60, fem_state_est_EKF(:,3),  'LineWidth', 2); hold on;
plot(time_vec/60, fem_states_IC_prop(:,3), 'LineWidth', 2);
plot(time_vec/60, fem_states(:,3),         'LineWidth', 2);
ylabel('$r_{z}\:\left(km\right)$','Interpreter','LaTex');
xlabel('Time (minutes)');
set(gca,'FontSize', 16, 'XLim',[0 max(time_vec/60)]);

% Plot 7: True positions of mothersat and femtosat
figure(7);
plot3(fem_states(:,1), fem_states(:,2), fem_states(:,3), 'linewidth', 4);
hold on;
plot3(mom_states(:,1), mom_states(:,2), mom_states(:,3), 'linewidth', 4);
[th, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(-pi, pi, 50));
[x,y,z]   = sph2cart(th, phi, rad_Titan*ones(50,50));
surface(x,y,z,'FaceColor', 'yellow');
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
axis equal;
set(gca,'FontSize', 16);

% Plot 8: position estimation error as a function of the angle between the
% femtosat and mothersat position vectors
figure(8);
for index = 1:m
    r = fem_states(index,1:3);
    R = mom_states(index,1:3);
    rR_angle(index) = subspace(r',R') * 180/pi;
end
plot(rR_angle, resids_pos_EKF,     'LineWidth', 2); hold on;
plot(rR_angle, resids_pos_IC_prop, 'LineWidth', 2);
% text(rR_angle(1)+2, resids_pos_IC_prop(1), 'Start', 'FontSize', 16);
% text(rR_angle(end), resids_pos_IC_prop(end)-5, 'End',   'FontSize', 16);
xlabel('Angle Between r and R (degrees)');
ylabel('Position Estimate Error (km)');
legend('EKF Estimate', 'IC Solution Propagation', 'Location', 'Best');
set(gca,'FontSize', 16);

% End timer
t_end = toc(t_start);
fprintf('Simulation run time: %d minutes and %i seconds\n', ...
        floor(t_end/60), floor(rem(t_end,60)));
