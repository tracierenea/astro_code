function [x_hat, fem_states] = iterated_EKF_orbit_det(time, x0_guess, y_meas, R, max_count, mom_states, freq, mu)
%ITERATED_EKF_ORBIT_DET Use extended Kalman filter to find IC
%   This is an implementation of the algorithm described in:
%   Optimal Estimation of Dynamic Systems, 2nd Edition
%   by Crassidis & Junkins
%   As described on page 476: "This procedure uses the extended
%   Kalman filter shown in Table 3.9 with Q = 0 to process the data
%   forward with some initial condition guess, and then process the
%   data backward to epoch."
%
%   Of course, the method has been adapted for the particular
%   measurement data type under study in this simulation (i.e.
%   Doppler shift).
%
max_step = (time(end)-time(1))/20;        % [1]
step_0   = max_step/3;                    % [1]
odeopt   = odeset ('RelTol', 1e-12, ...
                   'AbsTol', 1e-12, ...
                   'InitialStep', step_0, ...
                   'MaxStep', max_step);

% Initialize parameters
n     = length(x0_guess);                 % Number of states
m     = size(y_meas, 1);                  % number of measurements
P0    = 1*eye(n);                         % Initial covariance matrix
X0    = [x0_guess;
         reshape(P0, numel(P0), 1)];      % Stack state and covar.
count = 1;

% Allocate memory for storage containers
state_estimates_forwards      = zeros(m, n);
state_estimates_backwards     = zeros(m, n);
P_estimates_forwards          = zeros(m, n);
P_estimates_backwards         = zeros(m, n);
state_estimates_forwards(1,:) = x0_guess';
P_estimates_forwards(1,:)     = diag(P0)';

fprintf('Entering main EKF loop, standby... ')

% Main loop: don't exceed maximum number of iterations
while count <= max_count

  fprintf ('%i/%i ', count, max_count);

  % Process forwards
  for k = 1:m-1

    % Propagate the state and covariance together
    time_span    = [time(k) time(k+1)];
    [~, results] = ode45(@EKF_sys_eqns, time_span, X0, odeopt, mu,n);

    % Estimated satellite state at time k+1
    x_est = results(end, 1:n)';

    % Covariance matrix at time k+1
    P_est = reshape(results(end, n+1:end), n, n);

    % Calculate the H matrix
    mom_state = mom_states(k,:)';
    H = get_H_matrix_EKF(x_est, mom_state, freq);

    % Calculate the gain matrix K
    K = P_est*H'*inv(H*P_est*H' + R);

    % Update the state and covariance matrix
    h = get_measurements(x_est', mom_state', freq);
    x_est = x_est + K*(y_meas(k+1,1) - h);
    P_est = (eye(n) - K*H)*P_est;

    % Form the initial condition for the next sample interval
    X0 = [x_est; reshape(P_est, numel(P_est), 1)];

    % Store values for plotting later
    state_estimates_forwards(k+1,:) = x_est';
    P_estimates_forwards(k+1,:) = diag(P_est)';

  end % forward process 

  % Final state from forward pass becomes initial condition for
  % backward pass
  X_f = x_est;
  state_estimates_backwards(m, :) = X_f';

  % Reset covariance & create initial X0 for backwards pass
  P_estimates_backwards(m,:) = diag(P0)';
  X0 = [X_f; reshape(P0, numel(P0), 1)]; 

  % Process backwards
  for k = m:-1:2

    % Propagate the state and covariance together
    time_span    = [time(k) time(k-1)];
    [~, results] = ode45(@EKF_sys_eqns, time_span, X0, odeopt, mu,n);

    % Estimated satellite state at time k+1
    x_est = results(end, 1:n)';

    % Covariance matrix at time k+1
    P_est = reshape(results(end, n+1:end), n, n);

    % Calculate the H matrix
    mom_state = mom_states(k,:)';
    H = get_H_matrix_EKF(x_est, mom_state, freq);

    % Calculate the gain matrix K
    K = P_est*H'*inv(H*P_est*H' + R);

    % Update the state and covariance matrix
    h = get_measurements(x_est', mom_state', freq);
    x_est = x_est + K*(y_meas(k,1) - h);
    P_est = (eye(n) - K*H)*P_est;

    % Form the initial condition for the next sample interval
    X0 = [x_est; reshape(P_est, numel(P_est), 1)];

    % Store values for plotting later
    state_estimates_backwards(k,:) = x_est';
    P_estimates_backwards(k,:)     = diag(P_est)';

  end % backwards process

  % Reset covariance and use estimate from this iteration on the next
  X0 = [x_est; reshape(P0, numel(P0), 1)];


  % Increment counter
  count = count + 1;

  % Clear these containers before looping again
  if count <= max_count
    clear state_estimates_forwards state_estimates_backwards;
    state_estimates_forwards      = zeros(m, n);
    state_estimates_backwards     = zeros(m, n);
    P_estimates_forwards          = zeros(m, n);
    P_estimates_backwards         = zeros(m, n);
    state_estimates_forwards(1,:) = x_est';
    P_estimates_forwards(1,:)     = diag(P0)';
  end

  % TODO: add a test for convergence, like if r_est is changing by less than 10 m and r_dot is changing by less than .1 km/sec, and if so, break. This is just an optimization thing - do at the end.

end % main loop

% The state_estimates_backwards first state will be zeros, because k
% is only iterated to 2. So, set that point (t=0) equal to the first state in state_estimates_forwards
state_estimates_backwards(1,:) = state_estimates_forwards(1,:);

% Return results
x_hat = x_est;
fem_states = state_estimates_backwards;

% Footnote [1]
% "InitialStep and MaxStep in ode23..ode78 are not computed (as like
% the solvers in Matlab do) but must always be given by the user. If
% no value is given (eg. your first trial) then simply (tf-t0)/10 is
% taken - this might be completely wrong and depends on which ODE
% system should be solved. 
