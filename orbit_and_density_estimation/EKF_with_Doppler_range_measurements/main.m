%  Main simulalation script: orbit and atmospheric density estimation
%  from Doppler shift and range measurements.
%
%  Author: Tracie Perez // UTA MAE PhD Student
%
%  The purpose of this simulation is to estimate the initial state
%  of the femtosatellite's orbit (r, r_dot) using Doppler shift
%  measurements made on the mothersat.
%
%  Note: These scripts are configured for MATLAB. Scripts in other
%  directories of this repo are configured for octave. Why? Well MATLAB
%  can complete this script in under a minute while octave takes 7 mintues.
%  Also the octave GUI was getting really glitchy and I needed a GUI for 
%  rapid development. See the README for how to configure this for octave.
%
%

close all; clear all; clc;

%%% Start timer
t_start    = tic;

%%% Configure the simulation by defining the following constants:
n          = 8;                   % Number of states
G          = 6.6742e-11 / 1000^3; % km^3/kg*s^2, Univ. grav. constant
m_Titan    = 1.3452e23;           % kg, mass of Titan
mu         = G*m_Titan;           % km^3/s^2, Titan's grav. parameter
rad_Titan  = 2575.5;              % km, radius of Titan
freq       = 5.65e9;              % 5.650 GHz, freq of transmitted signal
rad_mom    = 1500 + rad_Titan;    % km
rad_fem    = 200  + rad_Titan;    % km
deploy_v   = 2/1000;              % km/sec (deployment velocity)
t0         = 0;                   % initial time
dt         = 1;                   % seconds, time interval
m          = 4801;                % make this many measurements; full decay
mass_fem   = 0.1;                 % kg, mass of femtosat
C_D        = 2.0;                 % dimensionless, C_D for sphere
radius_fem = 0.01;                % m, radius of femtosat sphere; 1cm
radius_fem = radius_fem/1000;     % convert to units of km
A          = pi*(radius_fem^2);   % cross-sectional area
BC         = mass_fem/(C_D*A);    % ballistic coefficient  
a_true     = 5.2680;              % This coefficient is later scaled by 1e9
b_true     =-0.0523;              % As a STATE, this is later scalled by 100

%%% Tune the EKF by experimenting with the following parameters:

% Error on the initial condition of the state provided to the EKF
error_vec  = [ones(3,1); 1e-1*ones(3,1); 1e-2; 1e-3];
% error_vec  = zeros(8,1);

% Measurement covariance (high = tell the EKF to trust the measurements)
R          = eye(2);%*20; 

% q = Q(t), variance of process noise for the EKF
q_EKF      = [1e-8*ones(3,1); 10; 10; 10; 1e-6; 1e-6];

% q = Q(t) for the dynamics, i.e. two_body function
q_dynamics = 1e-4;

% Initial covariance matrix
% P0         = diag([1e3*ones(1,3) 1e2*ones(1,3) 1e4 1e3]);
P0 = 100*eye(8);

% Specify standard deviation of noise added to each measurement type
noise_std_Doppler = 1;            % Std dev of Doppler shift meas. noise
noise_std_range   = 1;            % Std dev of range measurement noise

% Tried initially but failed
% P0         = diag(error_vec.^2); 
% R          = noise_std^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up conditions for simulation
y_dot0_fem = sqrt(mu/rad_fem);    % km/sec (circular orbit)
y_dot0_mom = sqrt(mu/rad_mom);    % km/sec (circular orbit)
X0_fem     = [rad_fem; 
              0;
              0;
              0;
              y_dot0_fem;
              deploy_v];
X0_mom     = [rad_mom;
              0; 
              0; 
              0; 
              y_dot0_mom; 
              0];

% True initial state
X0_true    = [X0_fem; a_true; b_true];

% Initial state "guess"
guess_error= error_vec.*randn(8,1);
X0_guess   = X0_true + guess_error;

% Generate truth data
tf              = t0 + m-1*dt;               % seconds, final time
time_vec        = [t0:dt:tf]';               % time vector for analysis
[~, mom_states] = ode45(@two_body_EOM, time_vec, X0_mom, [], mu);
[~, fem_states] = ode45(@two_body_EOM_with_drag, time_vec, X0_fem, [], ...
                        mu, BC, rad_Titan, a_true, b_true, q_dynamics);

% Generate measurement data
meas_data       = get_measurements(time_vec, fem_states, mom_states, freq);
time_meas       = meas_data(:,1);
y_true          = meas_data(:,2:3);

% Add noise to true measurements
meas_data(:,2)  = meas_data(:,2) + noise_std_Doppler*randn(size(y_true,1),1);
meas_data(:,3)  = meas_data(:,3) + noise_std_range * randn(size(y_true,1),1);

% Print out some info for confirmation
print_sim_description(m, meas_data, error_vec, dt, noise_std_Doppler, ...
    noise_std_range, guess_error, X0_mom, X0_fem, X0_guess);

% Filter measurements with EKF
[x_est, fem_state_est_EKF, P_EKF, y_hat_store] =  EKF(time_vec, X0_guess, meas_data, ...
R, mom_states, freq, mu, rad_Titan, BC, q_EKF, P0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print results
alt_final_est = norm([fem_state_est_EKF(end,1);
                      fem_state_est_EKF(end,2);
                      fem_state_est_EKF(end,3)]) - rad_Titan;
alt_final_true= norm(fem_states(end,1:3)) - rad_Titan;
a_estimate    = x_est(7);
b_estimate    = x_est(8)/100;
a_guess       = X0_guess(7);
b_guess       = X0_guess(8);
guess_error_a = a_true - a_guess;
guess_error_b = b_true - b_guess;
est_error_a   = a_true - a_estimate;
est_error_b   = b_true - b_estimate;
fprintf('Estimated final altitude at end of simulation: %.1f km\n', alt_final_est);
fprintf('True final altitude at end of simulation: %.1f km\n', alt_final_true);
fprintf('Results from EKF:\n')
fprintf('Coefficient A:\n');
fprintf('\tTruth\t   = %.6f\n', a_true);
fprintf('\tGuess\t   = %.6f\t(error = %.16f)\n', a_guess, guess_error_a);
fprintf('\tEstimate   = %.6f\t(error = %.16f)\n', a_estimate, est_error_a);
fprintf('Coefficient B:\n');
fprintf('\tTruth\t   = %.6f\n', b_true);
fprintf('\tGuess\t   = %.6f\t(error = %.6f)\n', b_guess, guess_error_b);
fprintf('\tEstimate   = %.6f\t(error = %.6f)\n', b_estimate, est_error_b);

% Create plots

% Plot 1: true satellite positions (x and y)
figure(1);
plot(mom_states(:,1), mom_states(:,2), 'r-', 'linewidth', 4); hold on;
img = imread('Titan_cropped.png');             % Load image of Titan
scale = 1.08;
x = [-rad_Titan*scale rad_Titan*scale];
y = [rad_Titan*scale -rad_Titan*scale];
image(x, y, img);                              % Plot the image
axis xy;
plot(fem_states(:,1), fem_states(:,2), 'k-', 'linewidth', 4);
xlabel('x (km)', 'FontSize', 16);
ylabel('y (km)', 'FontSize', 16);
title('True Satellite Positions');
axis equal;

% Plot 2: measurement data: truth and synthetically created meas's
figure(2);
subplot(211);
plot(time_vec/60, meas_data(:,2), 'o', time_vec/60, y_true(:,1), 'k-');
hold on;
plot(time_vec(1:m-1)/60, y_hat_store(:,1), 'r-');
ylabel('Doppler Shift (Hz)', 'FontSize', 16);
legend('Measurements', 'Truth', 'Estimates (y-hat)');
subplot(212);
plot(time_vec/60, meas_data(:,3), 'o', time_vec/60, y_true(:,2), 'k-');
hold on;
plot(time_vec(1:m-1)/60, y_hat_store(:,2), 'r-');
xlabel('Time (minutes)' , 'FontSize', 16);
ylabel('Range (km)', 'FontSize', 16);
legend('Measurements', 'Truth', 'Estimates (y-hat)');

% Estimated and true states
figure; % Components of position vector
subplot(311); 
plot(time_vec, fem_state_est_EKF(:,1), 'r-'); hold on;
plot(time_vec, fem_states(:,1), 'k-');
title('Position');
ylabel('r_x');
legend('Estimate', 'Truth');
subplot(312);
plot(time_vec, fem_state_est_EKF(:,2), 'r-'); hold on;
plot(time_vec, fem_states(:,2), 'k-');
ylabel('r_y');
subplot(313);
plot(time_vec, fem_state_est_EKF(:,3), 'r-'); hold on;
plot(time_vec, fem_states(:,3), 'k-');
ylabel('r_z');

figure; % Components of velocity vector
subplot(311); 
plot(time_vec, fem_state_est_EKF(:,4), 'r-'); hold on;
plot(time_vec, fem_states(:,4), 'k-');
title('Velocity');
ylabel('r_x dot');
legend('Estimate', 'Truth');
subplot(312);
plot(time_vec, fem_state_est_EKF(:,5), 'r-'); hold on;
plot(time_vec, fem_states(:,5), 'k-');
ylabel('r_y dot');
subplot(313);
plot(time_vec, fem_state_est_EKF(:,6), 'r-'); hold on;
plot(time_vec, fem_states(:,6), 'k-');
ylabel('r_z dot');

figure; % Estimates of coefficients a and b
subplot(211);
plot(time_vec, fem_state_est_EKF(:,7), 'r-'); hold on;
plot(time_vec, a_true*ones(size(time_vec)), 'k-');
ylabel('Coefficient a');
legend('Estimate', 'Truth');
subplot(212);
plot(time_vec, fem_state_est_EKF(:,8), 'r-'); hold on;
plot(time_vec, b_true*100*ones(size(time_vec)), 'k-');
ylabel('Coefficient b');
legend('Estimate', 'Truth');

% figure; % Components of acceleration due to drag vector
% for index = 1:m+1
%     r_vec_est = fem_state_est_EKF(index, 1:3);
%     r_vec_true = fem_states(index, 1:3);
%     v_vec_est = fem_state_est_EKF(index, 4:6);
%     v_vec_true = fem_states(index, 4:6);
%     a_est     = fem_state_est_EKF(index, 7);
%     b_est     = fem_state_est_EKF(index, 8)/1000;
%     v_mag_est     = norm(v_vec_est);
%     v_mag_true     = norm(v_vec_true);
%     altitude_est  = norm(r_vec_est) - rad_Titan;
%     altitude_true = norm(r_vec_true) - rad_Titan;
%     density_est   = a_est*exp(b_est*altitude_est);
%     density_true   = a_true*exp((b_true/1000)*altitude_true);
%     drag_est(index,1:3) = (-1/(2*BC))*density_est*v_mag_est*v_vec_est;
%     drag_true(index,1:3) = (-1/(2*BC))*density_true*v_mag_true*v_vec_true;
% end
% subplot(311);
% plot(time_vec, drag_est(:,1), 'r-'); hold on;
% plot(time_vec, drag_true(:,1), 'k-');
% ylabel('drag_x');
% legend('Estimate', 'Truth');
% subplot(312);
% plot(time_vec, drag_est(:,2), 'r-');hold on;
% plot(time_vec, drag_true(:,2), 'k-');
% ylabel('drag_y');
% subplot(313);
% plot(time_vec, drag_est(:,3), 'r-');hold on;
% plot(time_vec, drag_true(:,3), 'k-');
% ylabel('drag_z');

% Plot: satellite positions
figure; 
plot3(fem_states(:,1), fem_states(:,2), fem_states(:,3), 'g-');  % truth
hold on;
plot3(fem_state_est_EKF(:,1),  fem_state_est_EKF(:,2),  ...
      fem_state_est_EKF(:,3),  'r-');
plot3(mom_states(:,1), mom_states(:,2), mom_states(:,3), 'k-');
legend('Femtosat Truth', 'Femtosat Estimate', 'Mothersat');
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
[th, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(-pi, pi, 50));
[x,y,z]   = sph2cart(th, phi, rad_Titan*ones(50,50));
surface(x,y,z,'FaceColor', 'yellow'); % add sphere as placeholder for Titan

% Plot: position/velocity estimate error
fig = figure;
x_pos_resids = fem_state_est_EKF(:,1) - fem_states(:,1);
y_pos_resids = fem_state_est_EKF(:,2) - fem_states(:,2);
z_pos_resids = fem_state_est_EKF(:,3) - fem_states(:,3);
x_dot_resids = fem_state_est_EKF(:,4) - fem_states(:,4);
y_dot_resids = fem_state_est_EKF(:,5) - fem_states(:,5);
z_dot_resids = fem_state_est_EKF(:,6) - fem_states(:,6);
for index = 1:m
    resids_pos(index) = norm([x_pos_resids(index,1);
                              y_pos_resids(index,1);
                              z_pos_resids(index,1)]);
    resids_vel(index) = norm([x_dot_resids(index,1);
                              y_dot_resids(index,1);
                              z_dot_resids(index,1)]);
end
left_color = [0 0 153/255];
right_color = [0 102/255 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left;
plot(time_vec/60, resids_pos, '-', 'Linewidth', 2);
xlabel('Time (minutes)');
ylabel('Position Estimate Error (km)');
yyaxis right
plot(time_vec/60, resids_vel, '--', 'Linewidth', 2);
ylabel('Velocity Estimate Error (km/sec)');

% Combo plot  - % fix this for both measurement types.....
% figure;
% yyaxis left;
% plot(time_vec/60, meas_data(:,2), 'o', time_vec/60, y_true, 'k-');
% xlabel('Time (minutes)');
% ylabel('Doppler Shift (Hz)');
% legend('Measurements','Truth');
% yyaxis right
% plot(time_vec/60, resids_pos, 'r-');
% ylabel('Position Estimate Error (km)');

 
% % Plot 5: position error as a function of...geometric relationship, lat/lo
% % or maybe direction cosine angle?

% % Plot 6: femtosat altitude vs. time (i.e. decay)
figure;
counter = 1;
for index = 1:100:m
    r_true          = norm([fem_states(index,1),        ...
                            fem_states(index,2),        ...
                            fem_states(index,3)]);
    r_est           = norm([fem_state_est_EKF(index,1), ...
                            fem_state_est_EKF(index,2), ...
                            fem_state_est_EKF(index,3)]);
    z_true(counter) = r_true - rad_Titan;
    z_est(counter)  = r_est  - rad_Titan;
    counter = counter + 1;
end
plot(time_vec(1:100:m)/60, z_true, 'k-'); hold on;
plot(time_vec(1:100:m)/60, z_est,  'r-');
legend('Truth', 'Estimate');
xlabel('Time (minutes)');
ylabel('Altitude (km)');
title('Orbital Decay of Femtosatellite');

figure; % Estimated covariances
subplot(811); plot(P_EKF(:,1)); ylabel('P_{11} (r_x)');
subplot(812); plot(P_EKF(:,2)); ylabel('P_{22} (r_y)');
subplot(813); plot(P_EKF(:,3)); ylabel('P_{33} (r_z)');
subplot(814); plot(P_EKF(:,4)); ylabel('P_{44} (r_x dot)');
subplot(815); plot(P_EKF(:,5)); ylabel('P_{55} (r_y dot)');
subplot(816); plot(P_EKF(:,6)); ylabel('P_{66} (r_z dot)');
subplot(817); plot(P_EKF(:,7)); ylabel('P_{77} (a)');
subplot(818); plot(P_EKF(:,8)); ylabel('P_{88} (b)');
xlabel('Time (seconds)');

% End timer & print status
t_end = toc(t_start);
fprintf('Simulation run time: %d minutes and %i seconds\n', ...
        floor(t_end/60), floor(rem(t_end,60)));
