function [] = print_results(x0_true, x0_guess, x_estimate)
%PRINT_RESULTS Output final estimate of initial femtosat state
%   Inputs are:
%       x0_true    : the true initial state
%       x0_guess   : the initial guess provided to the estimator algorithm
%       x_estimate : the resulting initial state estimate from the
%                    NLS algorithm

% Residual = estimate - truth
error_vec = x_estimate - x0_true;

fprintf('Initial state of femtosat : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]  (truth)\n', x0_true); 

fprintf('Initial state guess       : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', x0_guess); 

fprintf('Final estimate of IC      : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', x_estimate);

fprintf('Estimate of IC error      : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]', error_vec);

fprintf('  %.1f km (pos error),',        norm(error_vec(1:3)));
fprintf('  %.3f km/sec (vel error)\n\n', norm(error_vec(4:6)));