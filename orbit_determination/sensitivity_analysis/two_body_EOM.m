function [Xdot] = two_body_EOM(~, X, mu)

  % Identify the state variables
  x      = X(1);
  y      = X(2);
  z      = X(3);
  x_dot  = X(4);
  y_dot  = X(5);
  z_dot  = X(6);

  % Apply 2-body orbital EOM
  r_mag    = norm([x y z]);
  x_ddot = -(mu/(r_mag^3))*x;
  y_ddot = -(mu/(r_mag^3))*y;
  z_ddot = -(mu/(r_mag^3))*z;

  % Return derivative of state vector
  Xdot = [x_dot; y_dot; z_dot; x_ddot; y_ddot; z_ddot];