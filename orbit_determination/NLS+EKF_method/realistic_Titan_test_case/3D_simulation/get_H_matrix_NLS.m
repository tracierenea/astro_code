function [H] = get_H_matrix_NLS(results, mom_states, freq)
%GET_H_MATRIX_NLS Return H matrix based on current estimate of IC
%   Evaluate the current estimate of the initial state of the
%   satellite by returning the corresponding H matrix.
%
%   Reference: Section 6.4, with the Psi matrix as 0 for now
%              Optimal Estimation of Dynamic Systems, 2nd Edition
%              by Crassidis & Junkins 

m     = size(results, 1);                    % number of measurements
n     = size(mom_states,2);                  % number of states
c     = 3e8;                                 % speed of light, m/s
lamda = (c/freq)/1000;                       % wavelength, km
n_phi = n^2;                                 % Phi is n x n

for ii = 1:m
  r            = results(ii,1:n/2)';
  r_dot        = results(ii,n/2+1:n)';
  R            = mom_states(ii,1:n/2)';
  R_dot        = mom_states(ii,n/2+1:n)';
  rho          =  r - R;                    % Range between mother and fem
  rho_mag      =  norm(rho);
  rho_dot_mag  =  (1/rho_mag)*(r'*r_dot - r'*R_dot - R'*r_dot + R'*R_dot);
  rho          = r-R;                        % column vector
  rho_mag      = norm(rho);
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
  
  Phi          = reshape(results(ii,n+1:n+n_phi), n, n);
  
  H(ii,:)      = dh_dx*Phi;                        % Eqn. 6.56
end