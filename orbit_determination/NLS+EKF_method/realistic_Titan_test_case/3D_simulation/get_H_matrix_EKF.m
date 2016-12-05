function [H] = get_H_matrix_EKF(fem_state, mom_state, freq)
%GET_H_MATRIX_EKF Return H matrix based on current estimate of IC
%   Evaluate the current estimate of the initial state of the
%   satellite by returning the corresponding H matrix.
%
%   Reference: Section 6.4 (but with modified measurement type)
%              Optimal Estimation of Dynamic Systems, 2nd Edition
%              by Crassidis & Junkins
%

n            = size(mom_state,1);           % number of states
c            = 3e8;                         % speed of light, m/s
lamda        = (c/freq)/1000;               % wavelength, km
r            = fem_state(1:n/2,   1);
R            = mom_state(1:n/2,   1);
r_dot        = fem_state(n/2+1:n, 1);
R_dot        = mom_state(n/2+1:n, 1);
rho          =  r - R;                      % Range between mother and fem
rho_mag      =  norm(rho);
rho_dot_mag  =  (1/rho_mag)*(r'*r_dot - r'*R_dot - R'*r_dot + R'*R_dot);
inv_rhoM     = 1/rho_mag;
inv_rhoM3    = 1/(rho_mag^3);  
temp         = r'*r_dot - r'*R_dot - R'*r_dot + R'*R_dot;
drhod_drx    = inv_rhoM*(r_dot-R_dot)'*[1;0;0] - ...
               inv_rhoM3*temp*(rho'*[1;0;0]);
drhod_dry    = inv_rhoM*(r_dot-R_dot)'*[0;1;0] - ...
               inv_rhoM3*temp*(rho'*[0;1;0]);
drhod_drz    = inv_rhoM*(r_dot-R_dot)'*[0;0;1] - ...
               inv_rhoM3*temp*(rho'*[0;0;1]);        
drhod_drxdot = inv_rhoM*(rho'*[1;0;0]);
drhod_drydot = inv_rhoM*(rho'*[0;1;0]);
drhod_drzdot = inv_rhoM*(rho'*[0;0;1]);
dh_dx        = (-1/lamda) * ...
               [drhod_drx    drhod_dry    drhod_drz      ...
                drhod_drxdot drhod_drydot drhod_drzdot];
dh_dx        = (-1/lamda) * ...
               [drhod_drx    drhod_dry    drhod_drz      ...
                drhod_drxdot drhod_drydot drhod_drzdot];
H            = dh_dx;                     % Table 3.9