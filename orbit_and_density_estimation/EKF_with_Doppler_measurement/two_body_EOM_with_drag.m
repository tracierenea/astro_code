function [Xdot] = two_body_EOM_with_drag(t,X,mu,BC,rad_Titan,a_true,b_true,q)

% Identify the state variables
x          =  X(1);
y          =  X(2);
z          =  X(3);
x_dot      =  X(4);
y_dot      =  X(5);
z_dot      =  X(6);

% Define position and velocity vectors
r_vec      =  [x; y; z];
r_mag      =  norm(r_vec);
v_vec      =  [x_dot; y_dot; z_dot];
v_mag      =  norm(v_vec);

% Calculate acceleration due to drag. Assume v_relative = -v
altitude   = r_mag - rad_Titan;
density    = (a_true*1e+09)*exp(b_true*altitude);
a_drag     = -0.5*(1/BC)*density*v_mag*v_vec;

% Apply 2-body orbital EOM
r_ddot_vec = (-mu/(r_mag^3))*r_vec + a_drag + q*randn(3,1);

% Return derivative of state vector
x_ddot     = r_ddot_vec(1);
y_ddot     = r_ddot_vec(2);
z_ddot     = r_ddot_vec(3);
Xdot       = [x_dot; y_dot; z_dot; x_ddot; y_ddot; z_ddot];