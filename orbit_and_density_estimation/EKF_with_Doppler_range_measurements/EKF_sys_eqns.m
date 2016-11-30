function big_vec_dot = EKF_sys_eqns(~,big_vector,mu,n,rad_planet,BC,q)
%EKF_SYS_EQNS System equations for propagation step in EKF
%   This function returns the derivative of the state vector and the
%   P matrix (which together form the "big state vector."
%
%   Reference: Table 3.9, page 188
%              Optimal Estimation of Dynamic Systems, 2nd Edition
%              by Crassidis & Junkins
%
%   For now we assume that the relative velocity is equal to -1 x the
%   femtosat's velocity vector, i.e. drag acts in the direction opposite
%   the velocity vector

% Identify the states and the P matrix
r_x       =  big_vector(1);
r_y       =  big_vector(2);
r_z       =  big_vector(3);
r_x_dot   =  big_vector(4);
r_y_dot   =  big_vector(5);
r_z_dot   =  big_vector(6);
a         =  big_vector(7);
% b         =  big_vector(8);
P         =  reshape(big_vector(n+1:end,1),n,n);
a_scaled  =  a*1e9;
b         =  big_vector(8)/100;

% Define position and velocity vectors
r_vec     =  [r_x; r_y; r_z];
r_mag     =  norm(r_vec);
altitude  =  r_mag - rad_planet;
v_vec     =  [r_x_dot; r_y_dot; r_z_dot];
v_mag     =  norm(v_vec);

% Solve for components of the acceleration
density   = a_scaled*exp(b*altitude);
drag_acc  = (-1/(2*BC))*density*v_mag*v_vec;
a_vec     = (-mu/(r_mag^3))*r_vec + drag_acc;

% Form f vector:
f         = [v_vec;
             a_vec;
             0;             % a is modeled as a constant, derivative = 0
             0];            % b is modeled as a constant; derviative = 0
         
G          = [0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0;
              0 0 0 0 0 0 0 0;
              0 0 0 1 0 0 0 0;
              0 0 0 0 1 0 0 0;
              0 0 0 0 0 1 0 0;
              0 0 0 0 0 0 1 0;
              0 0 0 0 0 0 0 1];
          
% Process noise 
Q         = diag(q);

% Form F matrix
F         = zeros(n,n);
F(1,4)    = 1;
F(2,5)    = 1;
F(3,6)    = 1;

temp      = (v_mag/(2*r_mag*BC))*a_scaled*exp(b*altitude);
F(4,1)    = 3*mu*(r_x^2)/(r_mag^5) - mu/(r_mag^3) - r_x*r_x_dot*temp;
F(5,2)    = 3*mu*(r_y^2)/(r_mag^5) - mu/(r_mag^3) - r_y*r_y_dot*temp;
F(6,3)    = 3*mu*(r_z^2)/(r_mag^5) - mu/(r_mag^3) - r_z*r_z_dot*temp;
F(4,2)    = 3*mu*r_x*r_y/(r_mag^5) - r_y*r_x_dot*temp; 
F(4,3)    = 3*mu*r_x*r_z/(r_mag^5) - r_z*r_x_dot*temp;
F(5,1)    = 3*mu*r_x*r_y/(r_mag^5) - r_x*r_y_dot*temp;
F(5,3)    = 3*mu*r_x*r_z/(r_mag^5) - r_z*r_y_dot*temp;
F(6,1)    = 3*mu*r_x*r_z/(r_mag^5) - r_x*r_z_dot*temp;
F(6,2)    = 3*mu*r_y*r_z/(r_mag^5) - r_y*r_z_dot*temp;

temp      = (-1/(2*BC))*a_scaled*exp(b*altitude);
F(4,4)    = temp*(v_mag + (r_x_dot^2)/v_mag);
F(5,5)    = temp*(v_mag + (r_y_dot^2)/v_mag);
F(6,6)    = temp*(v_mag + (r_z_dot^2)/v_mag);
F(4,5)    = temp*r_x_dot*r_y_dot*(1/v_mag);
F(4,6)    = temp*r_x_dot*r_z_dot*(1/v_mag);
F(5,4)    = temp*r_y_dot*r_x_dot*(1/v_mag);
F(5,6)    = temp*r_y_dot*r_z_dot*(1/v_mag);
F(6,4)    = temp*r_z_dot*r_x_dot*(1/v_mag);
F(6,5)    = temp*r_z_dot*r_y_dot*(1/v_mag);

temp      = (-1/(2*BC))*1e9*exp(b*altitude)*v_mag; %%%% YOOOOOOOOOOOOOOOO
% temp      = (-1/(2*BC))*exp(b*altitude)*v_mag;
F(4,7)    = temp*r_x_dot;
F(5,7)    = temp*r_y_dot;
F(6,7)    = temp*r_z_dot;

temp      = (-1/(2*BC))*a_scaled*v_mag*altitude*exp(altitude*b);
F(4,8)    = temp*r_x_dot;
F(5,8)    = temp*r_y_dot;
F(6,8)    = temp*r_z_dot;

% Propagation step
x_hat_dot = f;
P_dot     = F*P + P*F' + G*Q*G';

% Output derivative of the big state vector
big_vec_dot = [x_hat_dot;
               reshape(P_dot,n^2,1)];