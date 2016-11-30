function [x_hat, state_estimates, P_estimates, y_hat_store] = EKF(time, x0_guess, meas_data, R, mom_states, freq, mu, rad_planet, BC, q, P0)
%EKF Use extended Kalman filter to find atmospheric density model coeffic.
%   This is an implementation of the algorithm described in:
%
%   Optimal Estimation of Dynamic Systems, 2nd Edition
%   by Crassidis & Junkins
%   Section 7.3
%
%   Change made in implementation:
%   Rather than assume no process noise, q is provided as an input to this
%   function.
%
%   Also, the method has been adapted for the particular measurement data
%   type under study in this simulation (i.e. Doppler shift and range).
%
%   Inputs:
%     time            : array of time instances for mom states (i.e.
%                       the time interval to analyze)
%     x0_guess        : a guess for the initial state solution
%     meas_data       : 2-column matrix, [time   measurements]
%     R               : measurement covariance
%     max_count       : max # of iterations to allow
%     mom_states      : matrix of the mother satellite's states
%     freq            : transmission frequency
%     mu              : gravitational parameter of body being orbited
%     rad_planet      : radius of the planetary body being orbited
%     BC              : ballistic coefficient
%     q               : q = Q(t), variance of process noise
%     P0              : initial covariance matrix
%
%   Outputs:
%     x_hat           : estimate of initial femtosat state
%     state_estimates : matrix of femtosat state estimates
%     P_EKF           : error covariances
%

clear indices
time_meas= meas_data(:,1);  % Time for each measurement
y_meas_1 = meas_data(:,2);  % Doppler shift measurements
y_meas_2 = meas_data(:,3);  % Range measurements

% Scale b
x0_guess(8) = x0_guess(8)*100;

% Initialize parameters
n        = length(x0_guess);              % Number of states
m        = size(time_meas,1);             % Number of measurements
X0       = [x0_guess;
            reshape(P0, numel(P0), 1)];   % Stack state vector and P matrix

% Allocate memory for storage containers & initialize first entries
state_estimates      = zeros(m, n);
P_estimates          = zeros(m, n);
state_estimates(1,:) = x0_guess';
P_estimates(1,:)     = diag(P0)';
y_hat_store          = zeros(m-1,2);

% Loop over each time step in the interval
for index = 1:m-1

    % Propagate the state and covariance together
    time_span    = [time_meas(index) time_meas(index+1)];
    [~, results] = ode45(@EKF_sys_eqns,time_span,X0,[],mu,n,rad_planet,BC,q);

    % Estimated satellite state at time(index+1)
    x_est = results(end, 1:n)';

    % Covariance matrix at time(index+1)
    P_est = reshape(results(end, n+1:end), n, n);
    
    % The state of the mothersat is assumed known
    mom_state = mom_states(index+1,:)';

    % Calculate the H matrix at time(index+1)
    H = get_H_matrix_EKF(x_est, mom_state, freq);

    % Calculate the gain matrix K for time(index+1)
    K = P_est*H'/(H*P_est*H' + R);

    % Update the state and covariance matrix
    temp  = get_measurements(time(index+1), x_est(1:6)', mom_state', ...
                             freq); 
    h     = temp(2:3);
    y_k   = [y_meas_1(index,1); y_meas_2(index,1)];
    x_est = x_est + K*(y_k - h');
    P_est = (eye(n) - K*H)*P_est;
    
    % Store for debug later
    y_hat_store(index,:) = h;

    % Form the initial condition for the next sample interval
    X0 = [x_est; reshape(P_est, numel(P_est), 1)];

    % Store values for plotting later
    state_estimates(index+1,:) = x_est';
    P_estimates(index+1,:) = diag(P_est)';

end % for index = 1:m-1

% Return result
x_hat      = x_est;

