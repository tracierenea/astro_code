function [] = print_results(x0_true, x0_guess, x_estimate)
%PRINT_RESULTS Output final estimate of initial femtosat state
%   Inputs are:
%       x0_true    : the true initial state
%       x0_guess   : the initial guess provided to the estimator algorithm
%       x_estimate : the resulting initial state estimate from the
%                    NLS algorithm

error_vec = x0_true - x_estimate;

fprintf('Initial state of femtosat : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]  (truth)\n', x0_true); 

fprintf('Initial state guess       : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', x0_guess); 

fprintf('Final estimate of IC      : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', x_estimate);

fprintf('Estimate of IC error      : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n\n', error_vec);