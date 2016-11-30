function [x_hat, fem_states, P_EKF, y_hat_store] = EKF(time, x0_guess, meas_data, R, mom_states, freq, mu, rad_planet, BC, q, P0)
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
%   type under study in this simulation (i.e. Doppler shift).
%
%   Inputs:
%     time          : array of time instances for mom states (i.e.
%                     the time interval to analyze)
%     x0_guess      : a guess for the initial state solution
%     meas_data     : 2-column matrix, [time   measurements]
%     R             : measurement covariance
%     max_count     : max # of iterations to allow
%     mom_states    : matrix of the mother satellite's states
%     freq          : transmission frequency
%     mu            : gravitational parameter of body being orbited
%     rad_planet    : radius of the planetary body being orbited
%     BC            : ballistic coefficient
%     q             : q = Q(t), variance of process noise
%     P0            : initial covariance matrix
%
%   Outputs:
%     x_hat         : estimate of initial femtosat state
%     fem_states    : matrix of femtosat state estimates
%     P_EKF         : error covariances
%

clear indices
time_meas= meas_data(:,1);
y_meas   = meas_data(:,2);

% Initialize parameters
n        = length(x0_guess);              % Number of states
m        = size(y_meas, 1);               % Number of measurements
m_count  = 1;                             % Initialize a meas counters
% x0_guess(7) = x0_guess(7)*1e9;
x0_guess(8) = x0_guess(8)*100;
X0       = [x0_guess;
            reshape(P0, numel(P0), 1)];   % Stack state and P matrix
count    = 1;

% Identify the indices of the time vector when there are measurements
for index = 1:length(time_meas)
  indices(index,1) = find(time == time_meas(index));
end

% Allocate memory for storage containers
state_estimates      = zeros(m, n);
P_estimates          = zeros(m, n);
state_estimates(1,:) = x0_guess';
P_estimates(1,:)     = diag(P0)';
y_hat_store          = zeros(m-1,1);

% Process forwards; loop over each time step in the interval
for index = 1:size(time, 1)- 1

    % Propagate the state and covariance together
    time_span    = [time(index) time(index+1)];
    [~, results] = ode45(@EKF_sys_eqns,time_span,X0,[],mu,n,rad_planet,BC,q);

    % Estimated satellite state at time(index+1)
    x_est = results(end, 1:n)';

    % Covariance matrix at time(index+1)
    P_est = reshape(results(end, n+1:end), n, n);

    % Check if there is a measurement for this instant. If so,
    % invoke the state update step. If not, the state estimate will
    % be propagated forward
    if (m_count <= size(indices,1)) && (index == indices(m_count))

        % Calculate the H matrix at time(index+1)
        mom_state = mom_states(index+1,:)';
        H = get_H_matrix_EKF(x_est, mom_state, freq);

        % Calculate the gain matrix K for time(index+1)
        K = P_est*H'/(H*P_est*H' + R);

        % Update the state and covariance matrix
        temp  = get_measurements(time(index+1), x_est(1:6)', mom_state', ...
                                 freq); 
        h     = temp(2);
        y_hat_store(index) = h;
        x_est = x_est + K*(y_meas(m_count,1) - h);
        P_est = (eye(n) - K*H)*P_est;

        % increment counter
        m_count = m_count + 1;

    end % if there is a measurement for this instant

    % Form the initial condition for the next sample interval
    X0 = [x_est; reshape(P_est, numel(P_est), 1)];

    % Store values for plotting later
    state_estimates(index+1,:) = x_est';
    P_estimates(index+1,:) = diag(P_est)';

end % for index = 1:size(time, 1)- 1


% Return results
x_hat      = x_est;
fem_states = state_estimates;
P_EKF      = P_estimates;
