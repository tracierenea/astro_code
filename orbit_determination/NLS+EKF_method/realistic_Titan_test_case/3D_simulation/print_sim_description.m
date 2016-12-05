function [] = print_sim_description(test_case, y_meas, m_NLS, error_vec,...
    dt, noise_std, max_count, guess_error, X0_mom, X0_fem, x0_guess)
%PRINT_SIM_DESCRIPTION Print the configuration parameters for sim
%   Print out all of the parameters (either chosen or from random
%   number generators) the define the current simulation run

disp('')
fprintf('\tIterated extended Kalman filter for orbit determination');
fprintf(', test case: %i\n\n', test_case)
fprintf('First %i of the total %i measurements to be used for NLS\n', ...
        m_NLS, size(y_meas,1));
fprintf('Standard deviation of guess error  : ');
fprintf('[%.1f; %.1f; %.1f; %.2f; %.2f; %.2f]\n', error_vec)
fprintf('Time step (dt)                     : %.2f sec\n', dt);
fprintf('Std of noise specified to be       : %.2f Hz\n', noise_std);
fprintf('Number of iterations (NLS and EKF) : %i\n', max_count);
fprintf('Position vector guess error        : %.1f km\n', ...
        norm(guess_error(1:3)));
fprintf('Velocity vector guess error        : %.3f km/sec\n', ...
        norm(guess_error(4:6)));
fprintf('\n\t\t\t\tr_x      r_y     r_z')
fprintf('   r_x_dot r_y_dot r_z_dot (km or km/s)\n');
fprintf('Initial state of mothersat: ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]', X0_mom);
fprintf('  (truth)\n');
fprintf('Initial state of femtosat : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]', X0_fem);
fprintf('  (truth)\n');
fprintf('Initial state guess       : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', x0_guess);
fprintf('Initial state guess error : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', guess_error);