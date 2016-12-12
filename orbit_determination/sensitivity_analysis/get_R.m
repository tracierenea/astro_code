function R = get_R( t, n )
%GET_R Returns the R matrix as defined in Chapter 8 of Perez dissertation

% Get the full state transition matrix Phi
Phi  = get_STM( t, n );

% Grab each row
Phi1 = Phi(1,:);
Phi2 = Phi(2,:);
Phi3 = Phi(3,:);

% Compute R
R = Phi1'*Phi1 + Phi2'*Phi2 + Phi3'*Phi3;
end