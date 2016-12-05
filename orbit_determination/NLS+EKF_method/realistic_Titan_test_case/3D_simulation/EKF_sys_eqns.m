function big_vector_dot = EKF_sys_eqns(~, big_vector, mu, n)
%EKF_SYS_EQNS System equations for propagation step in EKF
%   This function returns the derivative of the state vector and the
%   P matrix by applying the 2-body orbital equation of motion.
%
%   Reference: Table 3.9, page 188 (with Q = 0)
%              Optimal Estimation of Dynamic Systems, 2nd Edition
%              by Crassidis & Junkins

  % Identify the state and the matrix vectors
  r_x      =  big_vector(1);
  r_y      =  big_vector(2);
  r_z      =  big_vector(3);
  r_x_dot  =  big_vector(4);
  r_y_dot  =  big_vector(5);
  r_z_dot  =  big_vector(6);
  P        =  reshape(big_vector(n+1:end,1),n,n);

  % Apply 2-body orbital EOM
  r_mag    =  norm([r_x r_y r_z]);
  r_x_ddot = -(mu/(r_mag^3))*r_x;
  r_y_ddot = -(mu/(r_mag^3))*r_y;
  r_z_ddot = -(mu/(r_mag^3))*r_z;

  % Calculate derivatives for Phi_dot
  F21_11   =  (3*mu*r_x^2)/(r_mag^5) - mu/r_mag^3;
  F21_22   =  (3*mu*r_y^2)/(r_mag^5) - mu/r_mag^3;
  F21_33   =  (3*mu*r_z^2)/(r_mag^5) - mu/r_mag^3;
  F21_12   =  (3*mu*r_x*r_y)/(r_mag^5);
  F21_13   =  (3*mu*r_x*r_z)/(r_mag^5);
  F21_21   =  (3*mu*r_x*r_y)/(r_mag^5);
  F21_23   =  (3*mu*r_y*r_z)/(r_mag^5);
  F21_31   =  (3*mu*r_x*r_z)/(r_mag^5);
  F21_32   =  (3*mu*r_y*r_z)/(r_mag^5);    
  F21      =  [F21_11   F21_12   F21_13;  ...
               F21_21   F21_22   F21_23;  ...    % Eqn. 6.58
               F21_31   F21_32   F21_33];
  F        =  [zeros(3) eye(3);  ...
               F21    zeros(3)];
  P_dot    =  F*P + P*F';                        % Table 3.9, Q = 0

  % Form the derivative of the big vector
  big_vector_dot = [r_x_dot;
                    r_y_dot;
                    r_z_dot;
                    r_x_ddot;
                    r_y_ddot;
                    r_z_ddot;
                    reshape(P_dot, numel(P_dot), 1)];