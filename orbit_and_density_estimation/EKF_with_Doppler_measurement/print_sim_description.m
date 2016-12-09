function [] = print_sim_description(m, y_meas, error_vec, dt, noise_std,...
                                    guess_error, X0_mom, X0_fem, x0_guess)
%PRINT_SIM_DESCRIPTION Print the configuration parameters for sim
%   Print out all of the parameters (either chosen or from random
%   number generators) the define the current simulation run

IC_pos_true   = X0_fem(1:3);
IC_vel_true   = X0_fem(4:6);
IC_pos_guess  = x0_guess(1:3);
IC_vel_guess  = x0_guess(4:6);

disp('')
fprintf('\tExtended Kalman filter for density model estimation\n\n');
fprintf('%i time points evaluated, %i measurements created', m, ...
        size(y_meas,1));
fprintf(' (%.1f %%)\n', (size(y_meas,1)/m)*100);
fprintf('Variance of guess error            : ');
fprintf('[%.1f; %.1f; %.1f; %.2f; %.2f; %.2f; %.2e; %.2e]\n', error_vec)
fprintf('Time step (dt)                     : %.2f  sec\n', dt);
fprintf('Simulation duration                : %.1f  min\n', m*dt/60);
fprintf('Std of noise specified to be       : %.3f Hz\n', noise_std);
fprintf('Position vector guess error        : %.3f km\n', ...
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
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', ...
        IC_pos_guess, IC_vel_guess);
fprintf('Initial state guess error : ');
fprintf('[%7.1f %7.1f %7.1f %7.3f %7.3f %7.3f]\n', ...
        x0_guess(1:6)-X0_fem);