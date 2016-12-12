function Phi = get_STM( t, n )
%GET_STM Calculate the 4 Clohessy-Wiltshire (CW) matrices
%   Calculate the state transition matrix, consisting of the 4 CW
%   submatrices, each of which is 3 x 3. These are defined in:
%       Orbital Mechanics for Engineering Students
%       Howard D. Curtis
%       Third Edition
%       Equations (7.53)
%
%   Inputs:
%   t   - time instant, in seconds
%   n   - mean motion of mothersat, in rad/sec
%
%   Output:
%   Phi - state transition matrix, 6 x 6
%

nt = n*t;

Phi_rr = [4-3*cos(nt)      0   0;                              % (7.53a)
          6*(sin(nt)-nt)   1   0;
          0                0   cos(nt)];
      
Phi_rv = [(1/n)*sin(nt)     (2/n)*(1-cos(nt))      0;          % (7.53b)
          (2/n)*(cos(nt)-1) (1/n)*(4*sin(nt)-3*nt) 0;
          0                  0                     (1/n)*sin(nt)];
          
Phi_vr = [3*n*sin(nt)     0   0;                               % (7.53c)
          6*n*(cos(nt)-1) 0   0;
          0               0  -n*sin(nt)];
      
Phi_vv = [ cos(nt)    2*sin(nt)    0;                          % (7.53d)
          -2*sin(nt)  4*cos(nt)-3  0;
           0          0            cos(nt)];

Phi = [Phi_rr   Phi_rv;                                        % (7.52)
       Phi_vr   Phi_vv];

end

