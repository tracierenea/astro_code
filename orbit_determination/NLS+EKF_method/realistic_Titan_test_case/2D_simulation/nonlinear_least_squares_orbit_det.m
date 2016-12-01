function [x_hat, fem_states] = nonlinear_least_squares_orbit_det(time, x0_guess, y_meas, Phi0_vec, invR, max_count, mom_states, freq, mu, radius_planet)

x_hat       = x0_guess;          % "hat" denotes current estimate
n           = length(x_hat);     % Number of states

% Initialize parameters just to get into the loop
count       = 1;
cost        = 1e6;     
last_cost   = 1e9;
cost_change = abs((cost - last_cost)/last_cost);

clear fem_states;

while (count <= max_count) && (cost_change > 0.01)

  % Integrate the system of equations to get states, Phi, and Psi
  big_vector       = [x_hat; Phi0_vec];
  [t_out, results] = ode45(@NLS_sys_eqns, time, big_vector, mu, n);
  fem_states       = results(:,1:4);

  % Obtain the predicted, or estimated, measurements based on the
  % current estimate
  temp  = get_measurements(time, fem_states, mom_states, freq, ...
                           radius_planet, 0);
  y_hat = temp(:,2);

  % Obtain the H matrix based on current estimate
  H = get_H_matrix_NLS(results, mom_states, freq);

  % Form the residuals
  resids = [y_meas - y_hat];

  % Compute updated cost
  cost = norm(resids);

  % Implement the update step
  del_x = inv(H'*invR*H)*H'*invR*resids;

  % Update the estimate
  x_hat = x_hat + del_x;

  % Increment the counter
  count = count + 1;
  cost_change = abs((cost - last_cost)/last_cost);
  last_cost = cost;

end