function [ X_dot ] = prop_function( ~, X, n, mu, rad_Titan, BC, q )
%PROP_FUNCTION Function for the propagation step in the main.m
%   This function is used in combination with ode45 to accomplish the 
%   propagate step of the continuous-discrete extended Kalman filter. 
%   Reference Table 3.9 on page 188 in Crassidis & Junkins.
%
%   Input X is the "big" state of the state vector x and the covariance
%   matrix P.
%

% Break apart the "big" state matrix
X_hat      = X(1:n);
P_vec      = X(n+1:end);

% Reshape the P_vec into the P_matrix
P          = reshape(P_vec,n,n);

% Break out the current estimate of each state in the state vector
r_x        = X_hat(1);
r_y        = X_hat(2);
r_z        = X_hat(3);
r_x_dot    = X_hat(4);
r_y_dot    = X_hat(5);
r_z_dot    = X_hat(6);
a          = X_hat(7);
b          = X_hat(8);
a_scaled   = a*1e9;
r_vec      = [r_x; r_y; r_z];
r_mag      = norm(r_vec);
altitude   = r_mag - rad_Titan;

% For now we assume that the relative velocity is equal to the velocity
% vector, i.e. drag acts in the direction opposite velocity
r_dot_vec  = [r_x_dot; r_y_dot; r_z_dot];
r_dot_mag  = norm(r_dot_vec);

% Solve for components of the acceleration
rho        = a_scaled*exp(b*altitude);
drag_accel = (-1/(2*BC))*rho*r_dot_mag*r_dot_vec;
r_ddot_vec = (-mu/(r_mag^3))*r_vec + drag_accel;

% Form f vector:
f          = [r_dot_vec;
              r_ddot_vec;
              0;             % a is modeled as a constant, derivative = 0
              0];            % b is modeled as a constant; derviative = 0
 
G          = [zeros(6,2);
              1 0
              0 1];
% Process noise 
Q          = q*eye(2);

% Form F matrix
F          = zeros(n,n);
F(1,4)     = 1;
F(2,5)     = 1;
F(3,6)     = 1;

temp       = (r_dot_mag/(2*r_mag*BC))*a_scaled*exp(b*altitude);
F(4,1)     = 3*mu*(r_x^2)/(r_mag^5) - mu/(r_mag^3) - r_x*r_x_dot*temp;
F(5,2)     = 3*mu*(r_y^2)/(r_mag^5) - mu/(r_mag^3) - r_y*r_y_dot*temp;
F(6,3)     = 3*mu*(r_z^2)/(r_mag^5) - mu/(r_mag^3) - r_z*r_z_dot*temp;
F(4,2)     = 3*mu*r_x*r_y/(r_mag^5) - r_y*r_x_dot*temp; 
F(4,3)     = 3*mu*r_x*r_z/(r_mag^5) - r_z*r_x_dot*temp;
F(5,1)     = 3*mu*r_x*r_y/(r_mag^5) - r_x*r_y_dot*temp;
F(5,3)     = 3*mu*r_x*r_z/(r_mag^5) - r_z*r_y_dot*temp;
F(6,1)     = 3*mu*r_x*r_z/(r_mag^5) - r_x*r_z_dot*temp;
F(6,2)     = 3*mu*r_y*r_z/(r_mag^5) - r_y*r_z_dot*temp;

temp       = (-1/(2*BC))*a_scaled*exp(b*altitude);
F(4,4)     = temp*(r_dot_mag + (r_x_dot^2)/r_dot_mag);
F(5,5)     = temp*(r_dot_mag + (r_y_dot^2)/r_dot_mag);
F(6,6)     = temp*(r_dot_mag + (r_z_dot^2)/r_dot_mag);
F(4,5)     = temp*r_x_dot*r_y_dot*(1/r_dot_mag);
F(4,6)     = temp*r_x_dot*r_z_dot*(1/r_dot_mag);
F(5,4)     = temp*r_y_dot*r_x_dot*(1/r_dot_mag);
F(5,6)     = temp*r_y_dot*r_z_dot*(1/r_dot_mag);
F(6,4)     = temp*r_z_dot*r_x_dot*(1/r_dot_mag);
F(6,5)     = temp*r_z_dot*r_y_dot*(1/r_dot_mag);

temp       = (-1/(2*BC))*1e9*exp(b*altitude)*r_dot_mag;
F(4,7)     = temp*r_x_dot;
F(5,7)     = temp*r_y_dot;
F(6,7)     = temp*r_z_dot;

temp       = (-1/(2*BC))*a_scaled*r_dot_mag*altitude*exp(altitude*b);
F(4,8)     = temp*r_x_dot;
F(5,8)     = temp*r_y_dot;
F(6,8)     = temp*r_z_dot;

% Propagation step
x_hat_dot  = f;
P_dot      = F*P + P*F' + G*Q*G';

% Output derivative of the big state vector
X_dot      = [x_hat_dot;
              reshape(P_dot,n^2,1)];

end


