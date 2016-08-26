function [x_hat, fem_states, P_f, P_b] = iterated_EKF_orbit_det(time, x0_guess, meas_data, R, max_count, mom_states, freq, mu, rad_planet, error_vec)
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
%   Inputs:
%     time          : array of time instances for mom states (i.e.
%                     the time interval to analyze)
%     x0_guess      : a guess for the initial state solution
%     meas_data     : 2-column matrix, [time   measurements]
%     R             : measurement covariance
%     max_count     : max # of iterations to allow
%     mom_states    : matrix of the mothersatellite's states
%     freq          : transmission frequency
%     mu            : gravitational parameter of body being orbitted
%     rad_planet    : radius of the planetary body being orbitted
%     error_vec     : used to define initial state guess in main.m
%
%   Output:
%     y             : 2-column matrix, [time measurements]
%

clear indices
time_meas= meas_data(:,1);
y_meas   = meas_data(:,2);
max_step = (time(end)-time(1))/20;        % see footnote [1]
step_0   = max_step/3;                    % see footnote [1]
odeopt   = odeset ('RelTol', 1e-12, ...
                   'AbsTol', 1e-12, ...
                   'InitialStep', step_0, ...
                   'MaxStep', max_step);

% Initialize parameters
n        = length(x0_guess);              % Number of states
m        = size(y_meas, 1);               % Number of measurements
m_count  = 1;                             % Initialize a meas counter
P0       = diag(error_vec.^2);            % Initial covariance matrix
X0       = [x0_guess;
            reshape(P0, numel(P0), 1)];   % Stack state and covar.
count    = 1;

% Identify the indices of the time vector when there are measurements
for index = 1:length(time_meas)
  indices(index,1) = find(time == time_meas(index));
end

% Allocate memory for storage containers
state_estimates_forwards      = zeros(m, n);
state_estimates_backwards     = zeros(m, n);
P_estimates_forwards          = zeros(m, n);
P_estimates_backwards         = zeros(m, n);
state_estimates_forwards(1,:) = x0_guess';
P_estimates_forwards(1,:)     = diag(P0)';

% Main loop: don't exceed maximum number of iterations
while count <= max_count

  % Process forwards; loop over each time step in the interval
  for index = 1:size(time, 1)- 1

    % Propagate the state and covariance together
    time_span    = [time(index) time(index+1)];
    [~, results] = ode45(@EKF_sys_eqns, time_span, X0, odeopt, mu,n);

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
      K = P_est*H'*inv(H*P_est*H' + R);

      % Update the state and covariance matrix. See footnote [2].
      h = get_measurements(time(index+1), x_est', mom_state', ...
                           freq, rad_planet, 0); 
      x_est = x_est + K*(y_meas(m_count,1) - h(2));
      P_est = (eye(n) - K*H)*P_est;

      % increment counter
      m_count = m_count + 1;

    end % if there is a measurement for this instant

    % Form the initial condition for the next sample interval
    X0 = [x_est; reshape(P_est, numel(P_est), 1)];

    % Store values for plotting later
    state_estimates_forwards(index+1,:) = x_est';
    P_estimates_forwards(index+1,:) = diag(P_est)';

  end % for index = 1:size(time, 1)- 1 // forward process

  % Final state from forward pass becomes initial condition for
  % backward pass
  X_f = x_est;
  state_estimates_backwards(m, :) = X_f';

  % Reset covariance & create initial X0 for backwards pass
  P_estimates_backwards(m,:) = diag(P0)';
  X0 = [X_f; reshape(P0, numel(P0), 1)]; 

  % Decrement counter
  m_count = m_count - 1;

  % Process backwards; loop over each time step in the interval
  for index = size(time, 1):-1:2
  % for index = m:-1:2

    % Propagate the state and covariance together
    time_span    = [time(index) time(index-1)];
    [~, results] = ode45(@EKF_sys_eqns, time_span, X0, odeopt, mu,n);

    % Estimated satellite state at time(index-1)
    x_est = results(end, 1:n)';
    
    % Covariance matrix at time(index-1)
    P_est = reshape(results(end, n+1:end), n, n);

    % Check if there is a measurement for this instant. If so,
    % invoke the state update step. If not, the state estimate will
    % be propagated forward
    if (m_count <= size(indices,1)) && (index == indices(m_count))

      % Calculate the H matrix
      mom_state = mom_states(index-1,:)';
      H = get_H_matrix_EKF(x_est, mom_state, freq);

      % Calculate the gain matrix K
      K = P_est*H'*inv(H*P_est*H' + R);

      % Update the state and covariance matrix. See footnote [3].
      h = get_measurements(time(index-1), x_est', mom_state', ...
                           freq, rad_planet, 0);
      x_est = x_est + K*(y_meas(m_count,1) - h(2));
      P_est = (eye(n) - K*H)*P_est;

      % increment counter
      m_count = m_count - 1;
    end % if there is a measurement for this instant

    % Form the initial condition for the next sample interval
    X0 = [x_est; reshape(P_est, numel(P_est), 1)];

    % Store values for plotting later
    state_estimates_backwards(index,:) = x_est';
    P_estimates_backwards(index,:)     = diag(P_est)';

  end % for index = size(time, 1):-1:2 // backwards process

  % Reset covariance because we have no new information. Use the
  % estimate from this iteration on the next one.
  X0 = [x_est; reshape(P0, numel(P0), 1)];

  % Increment counter
  count = count + 1;

  % Clear and reset these containers before looping again
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
x_hat      = x_est;
fem_states = state_estimates_backwards;
P_f        = P_estimates_forwards;
P_b        = P_estimates_backwards;


% Footnotes:
% [1]
% "InitialStep and MaxStep in ode23..ode78 are not computed (as like
% the solvers in Matlab do) but must always be given by the user. If
% no value is given (eg. your first trial) then simply (tf-t0)/10 is
% taken - this might be completely wrong and depends on which ODE
% system should be solved.
%
% [2]
% Note that this point, x_est and mom_state correspond to
% time(index+1). I incorrectly had time(index) previously.
%
% [3]
% Note that this point, x_est and mom_state correspond to
% time(index). I incorrectly had time(index+1) previously.
%
