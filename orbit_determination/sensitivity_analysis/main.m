%  Main script objective: Can the initial state be observed from Doppler
%  shift measurements?
%
%  Author: Tracie Perez // UTA MAE PhD Student
% 

close all; clear all; clc; format long;

%%% Start timer
t_start = tic;

% Because you can't pass additional parameters to fsolve..?
global n rho_dot_vec

%%% Configure the simulation by defining the following constants:
G          = 6.6742e-11 / 1000^3;  % km^3/kg*s^2, Univ. grav. constant
m_Titan    = 1.3452e23;            % kg, mass of Titan
mu         = G*m_Titan;            % km^3/s^2, Titan's grav. parameter
rad_Titan  = 2575.5;               % km, radius of Titan
freq       = 5.65e9;               % 5.650 GHz, freq of transmitted signal
rad_mom    = 1500+rad_Titan;       % km
n          = sqrt(mu/(rad_mom^3)); % rad/sec, mean motion of the mothersat
t          = [1:6];                % seconds, epochs of measurements
deploy_v   = 2/1000;               % km/sec (deployment velocity)
MC_runs    = 1000;                 % number of Monte Carlo runs to complete

% Standard deviation on the initial state guess error
error_vec = [10*ones(3,1); 0.01*ones(3,1)];

%%% Create true Doppler shift measurements, find true IC

%  1) specify true ICs

% mothersat on x-axis in circular orbit
y_dot0_mom = sqrt(mu/rad_mom);
x0_mom     = [rad_mom; 0; 0; 0; y_dot0_mom; 0];

% femtosat has just been deployed from the mothersat, so position is the
% same as mothersat, but velocity has z component
x0_fem     = x0_mom + [zeros(5,1); deploy_v];

% However, get states for one minute, and we're going to start our
% estimation process for the last 6 seconds of this interval (54 seconds
% after separation) so that Doppler isn't 0.

% 2) get true position and velocity states
time_vec        = 1:60;
[~, mom_states] = ode45(@two_body_EOM, time_vec, x0_mom, [], mu);
[~, fem_states] = ode45(@two_body_EOM, time_vec, x0_fem, [], mu);

% Discard all but the last 6 states
mom_states      = mom_states(55:end,:);     % 6x6
fem_states      = fem_states(55:end,:);     % 6x6

% These are now our true initial states
x0_mom          = mom_states(1,:)';
x0_fem          = fem_states(1,:)';

% 3) calculate rho_dot (not actually working with Doppler, but it's just
%    scaled by -1/lamda.
for index = 1:6
    r           = fem_states(index, 1:3)';  % femtosat position vector
    r_dot       = fem_states(index, 4:6)';  % femtosat velocity vector
    R           = mom_states(index, 1:3)';  % mothersat position vector
    R_dot       = mom_states(index, 4:6)';  % mothersat velocity vector
    rho         =  r - R;                   % range between mother and fem
    rho_mag     =  norm(rho);
    rho_dot_mag = (1/rho_mag)*(r'*r_dot-r'*R_dot-R'*r_dot+R'*R_dot);
    
    
    rho_dot_vec(index)   = rho_dot_mag;
    [r_rel, v_rel, ~]    = rva_relative(r, r_dot, R, R_dot, mu);
    rel_states(index,:)  = [r_rel' v_rel'];
    
end

% Solve for the true relative femtosat position and velocity in LVLH frame
r0_true = x0_fem(1:3); % true femtosat  position in inertial ref frame
v0_true = x0_fem(4:6); % true femtosat  velocity in inertial ref frame
R0_true = x0_mom(1:3); % true mothersat position in inertial ref frame
V0_true = x0_mom(4:6); % true mothersat velocity in inertial ref frame
[r0_true_rel, v0_true_rel, ~] = rva_relative(r0_true, v0_true, ...
                                             R0_true, V0_true, mu);

% Solve for the guess of the femtosat's position and velocity in the
%  LVLH frame
x0_true_rel  = [r0_true_rel; v0_true_rel];

% Solver was stopping prematurely at the default function evaluation limit
% of 600, so increase a bunch. Turn off diagnostics output for MC runs.
options = optimoptions('fsolve',                        ...
                       'MaxFunctionEvaluations', 10000, ...
                       'MaxIterations', 5000,           ...
                       'Diagnostics', 'off',            ...
                       'Display', 'off');
no_joy_counter = 0;
joy_counter    = 0;

for counter = 1:MC_runs
  fprintf('MC run: %i of %i\n', counter, MC_runs);
  
  guess_error  = error_vec.*randn(6,1);
  x0_guess_rel = x0_true_rel + guess_error;
  
  % Call fsolve to try to solve for the true relative position and velocity
  [x0_est_rel, fval, exitflag, output, jacobian] = fsolve(@fun, x0_guess_rel, options);
  
  if exitflag <= 0
      output.message
      no_joy_counter = no_joy_counter + 1;
      continue;
  else
      disp('fsolve status: equation solved.');
      joy_counter = joy_counter + 1;
  end
  
  % Save results
  guess_error = x0_true_rel - x0_guess_rel;
  est_error   = x0_true_rel - x0_est_rel;
  x0_guess_store(joy_counter,:)    = x0_guess_rel;
  x0_est_store(joy_counter,:)      = x0_est_rel;
  guess_error_store(joy_counter,:) = guess_error;
  est_error_store(joy_counter,:)   = est_error;
  
  % Clear variables for next iteration
  clear x0_guess_rel x0_est_rel;
  
end

figure; 
for index = 1:joy_counter
    plot3(guess_error_store(index,1), ...
          guess_error_store(index,2), ...
          guess_error_store(index,3), 'g^'); hold on;
    plot3(est_error_store(index,1),   ...
          est_error_store(index,2),   ...
          est_error_store(index,3), 'rv');
end
xlabel('$\delta x\:\left(km\right)$', 'Interpreter', 'LaTex');
ylabel('$\delta y\:\left(km\right)$', 'Interpreter', 'LaTex');
zlabel('$\delta z\:\left(km\right)$', 'Interpreter', 'LaTex');
set(gca,'FontSize',16);

% All of the error calculations below are in the LVLH frame
std_guesses = [std(guess_error_store(:,1)); ...
               std(guess_error_store(:,2)); ...
               std(guess_error_store(:,3)); ...
               std(guess_error_store(:,4)); ...
               std(guess_error_store(:,5)); ...
               std(guess_error_store(:,6))];
std_est     = [std(est_error_store(:,1));
               std(est_error_store(:,2));
               std(est_error_store(:,3));
               std(est_error_store(:,4));
               std(est_error_store(:,5));
               std(est_error_store(:,6))];
fprintf('\n\n\n');
fprintf('Std of pos guess errors    : %.3f  %.3f  %.3f\n', std_guesses(1:3));
fprintf('Std of pos estimate errors : %.3f  %.3f  %.3f\n', std_est(1:3)); 
fprintf('Std of vel guess errors    : %.3f  %.3f  %.3f\n', std_guesses(4:6));
fprintf('Std of vel estimate errors : %.3f  %.3f  %.3f\n', std_est(4:6));


% End timer
t_end = toc(t_start);
fprintf('\n\nScript run time: %d minutes and %i seconds\n', ...
        floor(t_end/60), floor(rem(t_end,60)));
    
% Save all variables from the current workspace
save main_sim.mat;