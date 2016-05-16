function [H] = get_H_matrix(results, mom_states, freq)
%GETHMATRIX Return H matrix based on current estimate of IC
%   Evaluate the current estimate of the initial state of the
%   satellite by returning the corresponding H matrix.

m     = size(results, 1);                    % number of measurements
n     = size(mom_states,2);                  % number of states
c     = 3e8;                                 % speed of light, m/s
lamda = (c/freq)/1000;                       % wavelength, km
n_phi = n^2;                                 % Phi is n x n

for ii = 1:m
  r            = [   results(ii,1);    results(ii,2)];
  R            = [mom_states(ii,1); mom_states(ii,2)];
  r_dot        = [   results(ii,3);    results(ii,4)];
  R_dot        = [mom_states(ii,3); mom_states(ii,4)];
  rho_dot_mag  = get_rho_dot_mag(r, r_dot, R, R_dot);
  rho          = r-R;
  rho_mag      = norm(rho);
  inv_rhoM     = 1/rho_mag;
  inv_rhoM3    = 1/(rho_mag^3);  
  
  temp         = r'*r_dot - r'*R_dot - R'*r_dot + R'*R_dot;
  
  drhod_drx    = inv_rhoM*(r_dot-R_dot)'*[1;0] - ...
                 inv_rhoM3*temp*(rho'*[1;0]);
  
  drhod_dry    = inv_rhoM*(r_dot-R_dot)'*[0;1] - ...
                 inv_rhoM3*temp*(rho'*[0;1]);
  
  drhod_drxdot = inv_rhoM*(rho'*[1;0]);
  
  drhod_drydot = inv_rhoM*(rho'*[0;1]);
  
  dh_dx        = (-1/lamda) * ...
                 [drhod_drx drhod_dry drhod_drxdot  drhod_drydot];
  
  Phi          = reshape(results(ii,n+1:n+n_phi), n, n);
  
  H(ii,:)      = dh_dx*Phi;                        % Eqn. 6.56
end