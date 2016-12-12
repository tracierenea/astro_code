function F = fun( x0 )
%FUN Nonlinear system of equations to solve
%   Nonlinear system of equations to solve to find the initial state of the
%   femtosat, that is, the relative position and velocity in the LVLH frame

global n rho_dot_vec

Q1 = get_Q(1,n);
Q2 = get_Q(2,n);
Q3 = get_Q(3,n);
Q4 = get_Q(4,n);
Q5 = get_Q(5,n);
Q6 = get_Q(6,n);

R1 = get_R(1,n);
R2 = get_R(2,n);
R3 = get_R(3,n);
R4 = get_R(4,n);
R5 = get_R(5,n);
R6 = get_R(6,n);

F(1) = (x0'*Q1*x0) / sqrt(x0'*R1*x0) - rho_dot_vec(1);
F(2) = (x0'*Q2*x0) / sqrt(x0'*R2*x0) - rho_dot_vec(2);
F(3) = (x0'*Q3*x0) / sqrt(x0'*R3*x0) - rho_dot_vec(3);
F(4) = (x0'*Q4*x0) / sqrt(x0'*R4*x0) - rho_dot_vec(4);
F(5) = (x0'*Q5*x0) / sqrt(x0'*R5*x0) - rho_dot_vec(5);
F(6) = (x0'*Q6*x0) / sqrt(x0'*R6*x0) - rho_dot_vec(6);

end

