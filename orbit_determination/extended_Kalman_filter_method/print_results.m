function [] = print_results(x0_true, x0_guess, x_estimate)
%PRINT_RESULTS Output final estimate of initial femtosat state
%   Inputs are:
%       x0_true    : the true initial state
%       x0_guess   : the initial guess provided to the NLS algorithm
%       x_estimate : the resulting initial state estimate from the
%                    NLS algorithm

fprintf('\n\t\t  r_x, km    r_y, km    r_x_dot   r_y_dot, km/s\n');

fprintf('True IC        :%9.1f  %9.1f  %9.4f  %9.4f\n',   ...
        x0_true(1), x0_true(2), x0_true(3), x0_true(4)); 

fprintf('Initial guess  :%9.1f  %9.1f  %9.4f  %9.4f\n',   ...
        x0_guess(1), x0_guess(2), x0_guess(3), x0_guess(4)); 

fprintf('Final estimate :%9.1f  %9.1f  %9.4f  %9.4f\n',   ...
        x_estimate(1), x_estimate(2), x_estimate(3), x_estimate(4));

fprintf('Error          :%9.1f  %9.1f  %9.4f  %9.4f\n\n', ...
        x0_true(1) - x_estimate(1), ...
        x0_true(2) - x_estimate(2), ...
        x0_true(3) - x_estimate(3), ...
        x0_true(4) - x_estimate(4));