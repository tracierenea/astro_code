function Q = get_Q( t, n )
%GET_Q Returns the Q matrix as defined in Chapter 8 of Perez dissertation

% Get the full state transition matrix Phi
Phi  = get_STM( t, n );

% Grab each row of the state transition matrix, Phi
Phi1 = Phi(1,:);
Phi2 = Phi(2,:);
Phi3 = Phi(3,:);
Phi4 = Phi(4,:);
Phi5 = Phi(5,:);
Phi6 = Phi(6,:);

% Compute Q as we've defined it
Q = Phi1'*Phi4 + Phi2'*Phi5 + Phi3'*Phi6;

end