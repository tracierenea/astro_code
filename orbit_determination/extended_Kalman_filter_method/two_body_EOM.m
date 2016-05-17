function [Xdot] = two_body_EOM(~, X, mu)

  % Identify the state and the matrix vectors
  x      = X(1);
  y      = X(2);
  x_dot  = X(3);
  y_dot  = X(4);

  % Apply 2-body orbital EOM
  r_mag    = norm([x y]);
  x_ddot = -(mu/(r_mag^3))*x;
  y_ddot = -(mu/(r_mag^3))*y;

  % Return derivative of state vector
  Xdot = [x_dot; y_dot; x_ddot; y_ddot];