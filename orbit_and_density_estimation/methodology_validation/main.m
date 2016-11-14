% Show that there is enough information in the position and velocity states
% to extract atmospheric density estimates. i.e. Apply an EKF to truth data
%
% Author: Tracie Perez // UTA MAE PhD Student
%

clc; close all; clear all;

% Constants for the simulation
n            = 8;                   % Number of states
G            = 6.6742e-11 / 1000^3; % km^3/kg*s^2, Univ. grav. constant
m_Titan      = 1.3452e23;           % kg, mass of Titan
mu           = G*m_Titan;           % km^3/s^2, Titan's grav. parameter
rad_Titan    = 2575.5;              % km, radius of Titan
t0           = 0;                   % initial time
dt           = 1;                   % seconds, time interval
m            = 4801;                % try to make this many meas's
mass_fem     = 0.1;                 % kg, mass of femtosat
C_D          = 2.0;                 % dimensionless, C_D for sphere
radius_fem   = 0.01;                % m, radius of femtosat sphere; 1cm
radius_fem   = radius_fem/1000;     % convert to km
A            = pi*(radius_fem^2);   % cross-sectional area
BC           = mass_fem/(C_D*A);    % ballistic coefficient
rad_fem      = 200+rad_Titan;       % km
y_dot0       = sqrt(mu/rad_fem);    % km/sec (circular orbit)
deploy_v     = 2/1000;              % km/sec (deployment velocity)
X0_fem       = [rad_fem; 
                0;
                0;
                0;
                y_dot0;
                deploy_v];
a_true       =  5.2680;             % coefficient scaled by e+09;
b_true       = -0.0523;

% Experimenting with values R and q
R            =  1*eye(6);           % R: the variance of measurement noise
q            =  1e-6;               % q = Q(t), variance of process noise

% Initial guess for P0
P0           = eye(n);
P0(7,7)      = 10;
P0(8,8)      = 10;

% Initial state guesses:
X0_fem_guess = X0_fem;
a_guess      = a_true*1.5;
b_guess      = b_true*1.5;

% The state is: X = [r_x r_y r_z r_x_dot r_y_dot r_z_dot a b]
X0_guess     = [X0_fem_guess; a_guess; b_guess];

% Create true state data
tf           = t0 + (m-1)*dt;       % seconds, final time
time_vec     = [t0:dt:tf]';         % time vector for analysis

% Create synthetic measurement data
[~,states]   = ode45(@two_body_EOM_with_drag, time_vec, X0_fem, [], mu, ...
                     BC, rad_Titan);

% Specify H matrix
H            = [1 0 0 0 0 0 0 0;
                0 1 0 0 0 0 0 0;
                0 0 1 0 0 0 0 0;
                0 0 0 1 0 0 0 0;
                0 0 0 0 1 0 0 0;
                0 0 0 0 0 1 0 0];

% Initialize variables
P_minus_store(:,:,1) = P0;
X_hat_minus(1,:) = X0_guess';

% Run the extended Kalman filter algorithm given in Table 3.9
for k = 1:m-1

    % Update gain
    P_minus         = P_minus_store(:,:,k);
    K_k             = P_minus*H'/(H*P_minus*H' + R);
    
    % Update step (for this problem, y_meas is the true state)
    X_hat_plus(k,:) = X_hat_minus(k,:)' + ...
                      K_k*(states(k,1:6)' - X_hat_minus(k,1:6)');
    P_plus          = (eye(n)-K_k*H)*P_minus;
    
    % Form the "big" state to propagate
    X0              = [X_hat_plus(k,:)'; 
                       reshape(P_plus,n^2,1)];
    
    % Propagate step
    [t_out,X_out]   = ode45(@prop_function, [time_vec(k) time_vec(k+1)],...
                            X0, [], n, mu, rad_Titan, BC, q);
    
    % Now the output becomes the new x_hat_minus and P_k_minus
    X_last                 = X_out(end,:);
    X_hat_minus(k+1,:)     = X_last(1:n);
    P_minus_vec            = X_last((n+1):end);
    P_minus                = reshape(P_minus_vec,n,n);
    P_minus_store(:,:,k+1) = P_minus;
    
    % Store the 3*sigma bounds
    sig3(k,:) = (3*diag(P_plus).^(0.5))';
    
end

% Output simulation results
alt_final     = norm([states(end,1), states(end,2), states(end,3)]) - ...
                rad_Titan;
error_A       = a_true - X_hat_minus(end,7);
error_B       = b_true - X_hat_minus(end,8);
guess_error_A = a_true - a_guess;
guess_error_B = b_true - b_guess;
A_estimate    = X_hat_minus(end,7);
B_estimate    = X_hat_minus(end,8);
est_error_A   = a_true - A_estimate;
est_error_B   = b_true - B_estimate;
fprintf('Initial state of femtosat: ');
fprintf('%.3f  %.3f  %.3f  %.3f  %.3f  %.3f\n',X0_fem)
fprintf('Duration of simulation: %i sec (%.1f min)\n', m*dt, m*dt/60);
fprintf('True final altitude at end of simulation: %.1f km\n', alt_final);
fprintf('Measurments processed: %i\n\n', k);
fprintf('Coefficient A:\n');
fprintf('\tTruth\t   = %.6f\n', a_true);
fprintf('\tGuess\t   = %.6f\t(error = %.6f)\n', a_guess, guess_error_A);
fprintf('\tEstimate   = %.6f\t(error = %.6f)\n', A_estimate, est_error_A);
fprintf('Coefficient B:\n');
fprintf('\tTruth\t   = %.6f\n', b_true);
fprintf('\tGuess\t   = %.6f\t(error = %.6f)\n', b_guess, guess_error_B);
fprintf('\tEstimate   = %.6f\t(error = %.6f)\n', B_estimate, est_error_B);


% Create plots
figure;
img = imread('Titan_cropped.png');             % Load image of Titan
scale = 1.08;
x = [-rad_Titan*scale rad_Titan*scale];
y = [rad_Titan*scale -rad_Titan*scale];
image(x, y, img);                                % Plot the image
axis xy;
hold on;
plot(states(:,1), states(:,2), 'k-', 'linewidth', 4);
axis equal;
xlabel('x (km)', 'FontSize', 16);
ylabel('y (km)', 'FontSize', 16);

% figure; 
% plot3(states(:,1), states(:,2), states(:,3), 'k-');     % truth
% xlabel('x (km)', 'FontSize', 14);
% ylabel('y (km)', 'FontSize', 14);
% zlabel('z (km)', 'FontSize', 14); 

figure;
for index = 1:m
    altitudes(index) = norm([states(index,1), states(index,2)])-rad_Titan;
end
plot(time_vec/60, altitudes, 'k-', 'linewidth', 3);
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Altitude (km)' , 'FontSize', 16)

figure;
a_estimate_resids = X_hat_minus(:,7) - a_true;
plot(time_vec(1:m-1)/60, -sig3(:,7), 'r-'); hold on;
plot(time_vec(1:m-1)/60,  sig3(:,7), 'r-');
plot(time_vec(1:m)  /60,  a_estimate_resids, 'k-');
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Coefficient ''a'' Residual', 'FontSize', 16);
% title('Test Case 1', 'FontSize', 16);

figure;
b_estimate_resids = X_hat_minus(:,8) - b_true;
plot(time_vec(1:m-1)/60, -sig3(:,8), 'r-'); hold on;
plot(time_vec(1:m-1)/60,  sig3(:,8), 'r-'); 
plot(time_vec(1:m)  /60,  b_estimate_resids, 'k-');
axis([-1 time_vec(end)/60 -0.12 0.12])
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Coefficient ''b'' Residual', 'FontSize', 16);
% title('Test Case 1', 'FontSize', 16);
 