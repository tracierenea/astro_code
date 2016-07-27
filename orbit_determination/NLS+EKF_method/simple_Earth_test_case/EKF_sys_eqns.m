function big_vector_dot = EKF_sys_eqns(~, big_vector, mu, n)
%EKF_SYSEQNS System equations for propagation step in EKF
%   This function returns the derivative of the state vector and the
%   P matrix by applying the 2-body orbital equation of motion.
%
%   Reference: Table 3.9, page 188 (with Q = 0)
%              Optimal Estimation of Dynamic Systems, 2nd Edition
%              by Crassidis & Junkins

  % Identify the state and the matrix vectors
  r_x      = big_vector(1);
  r_y      = big_vector(2);
  r_x_dot  = big_vector(3);
  r_y_dot  = big_vector(4);
  P        = reshape(big_vector(n+1:end,1),n,n);

  % Apply 2-body orbital EOM
  r_mag    = norm([r_x r_y]);
  r_x_ddot = -(mu/(r_mag^3))*r_x;
  r_y_ddot = -(mu/(r_mag^3))*r_y;

  % Calculate derivatives for Phi_dot
  F21_11   = (3*mu*r_x^2)/(r_mag^5) - mu/r_mag^3;
  F21_22   = (3*mu*r_y^2)/(r_mag^5) - mu/r_mag^3;
  F21_12   = (3*mu*r_x*r_y)/(r_mag^5);
  F21_21   = (3*mu*r_x*r_y)/(r_mag^5);
  F21      = [F21_11   F21_12;  ...
              F21_21   F21_22];              % Eqn. 6.58
  F        = [zeros(2) eye(2);  ...
              F21    zeros(2)];              % Eqn. 6.57
  P_dot    = F*P + P*F';                     % Table 3.9, Q = 0

  % Form the derivative of the big vector
  big_vector_dot = [r_x_dot;
                    r_y_dot;
                    r_x_ddot;
                    r_y_ddot;
                    reshape(P_dot, numel(P_dot), 1)];

