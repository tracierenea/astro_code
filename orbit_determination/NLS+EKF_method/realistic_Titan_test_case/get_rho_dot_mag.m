function [ rho_dot_mag ] = get_rho_dot_mag(r, r_dot, R, R_dot)
%GET_RHO_DOT_MAG Return magnitude of rho_dot vector
%   Inputs:
%   r     = position vector of femtosat
%   R     = position vector of mothersat
%   r_dot = velocity vector of femtosat
%   R_dot = velocity vector of mothersat

  rho         = r - R;
  rho_mag     = norm(rho);
  rho_dot_mag = (1/rho_mag)*(r'*r_dot - r'*R_dot - ...
                             R'*r_dot + R'*R_dot);
end

