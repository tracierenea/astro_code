%  Author: Tracie Perez
%
% The purpose of this script is to evaluate the Fisher information
% matrix, as a way to understand the information available in the 
% Doppler shift measurements used in this simulation.

% F = H'*inv(R)*H
%   = a lower bound on the expected errors between the estimated
%     quantities and the true values from the known properties of
%     the measurement errors.
%
close all; clear all; clc;

% Start timer
tic;

% Pick a case study
% test_case = 1;
test_case = 2;
% test_case = 3;

% If there are warnings, find out which functions generated them
debug_on_warning (true);

% This is the true mothersat initial state; it's the same in all 3
% test cases
G          = 6.6742e-11 / 1000**3;        % km^3/kg*s^2, uni g const 
m_Titan    = 1.3452e23;                   % kg, mass of Titan
mu         = G*m_Titan;                   % km^3/s^2, Titan's g param
rad_Titan  = 2575.5;                      % km, radius of Titan
freq       = 5.65e9;                      % 5.650 GHz
rad_mom    = 1500+rad_Titan;              % km
y_dot0_mom = sqrt(mu/rad_mom);            % km/sec (circular orbit)
mom_state0 = [rad_mom; 0; 0; y_dot0_mom]; % truth
deploy_v   = 2/1000;              % km/sec (deployment velocity)

% Standard deviation of noise, as specified in main.m
if     test_case == 1
  noise_std   = 10;
  rad_fem     = rad_mom;
  y_dot0_fem  = sqrt(mu/rad_fem);    % km/sec (circular orbit)
  X0_fem      = [rad_fem; 0; -deploy_v; y_dot0_fem]; % truth
  KF_MC_results = load('previousMCresults/MC_results_case1.txt');
elseif test_case == 2
  noise_std   = 500;
  rad_fem     = 400 + rad_Titan;     % km
  y_dot0_fem  = sqrt(mu/rad_fem);    % km/sec (circular orbit)
  X0_fem      = [rad_fem; 0; 0; y_dot0_fem];
  KF_MC_results = load('previousMCresults/MC_results_case2.txt');
elseif test_case == 3
  noise_std   = 10;
  rad_fem     = 400 + rad_Titan;     % km
  y_dot0_fem  = sqrt(mu/rad_fem);    % km/sec (circular orbit)
  vel_fem     = y_dot0_fem;
  vel_mom     = y_dot0_mom;
  X0_fem      = [ rad_fem*cosd(60);
                  rad_fem*sind(60);
                 -vel_fem*cosd(30);
                  vel_fem*sind(30)];
  X0_mom      = [ rad_mom*cosd(30);
                  rad_mom*sind(30);
                 -vel_mom*cosd(60);
                  vel_mom*sind(60)];
  KF_MC_results = load('previousMCresults/MC_results_case3.txt');
end

% These are the estimates from the Monte Carlo runs
estimates = KF_MC_results(:,10:13);

% Compute inverse of measurement covariance matrix just once
R        = noise_std^2;             % Measurement covariance
invR     = 1/R;

test_run_vec = [1:size(estimates,1)];
for MC_count = test_run_vec
  % Femtosat initial state estimate results for this test run
  this_estimate = estimates(MC_count,:)';

  % Obtain the H for this MC run 
  H = get_H_matrix_EKF(this_estimate, mom_state0, freq);

  % Calculate error
  error_vec(1:4,MC_count) = this_estimate - X0_fem;
  
  F = H'*invR*H;
%  P = inv(F);  % <---- problem not numerically conditioned to calc P
  sigma_x(MC_count) = F(1,1)^0.5;
  sigma_y(MC_count) = F(2,2)^0.5;
  sigma_xdot(MC_count) = F(3,3)^0.5;
  sigma_ydot(MC_count) = F(4,4)^0.5;

end
figure();
%plot(test_run_vec, error_vec(1,:),'*', test_run_vec, 3*sigma_x,'-',test_run_vec,-3*sigma_x,'-');
plot(test_run_vec, 3*sigma_x,'-');
figure();
%plot(test_run_vec, error_vec(2,:),'*', test_run_vec, 3*sigma_y,'-',test_run_vec,-3*sigma_y,'-');
plot(test_run_vec, 3*sigma_y,'-');
figure();
%plot(test_run_vec, error_vec(3,:),'*', test_run_vec, 3*sigma_xdot,'-',test_run_vec,-3*sigma_xdot,'-');
plot(test_run_vec, 3*sigma_xdot,'-');
figure();
%plot(test_run_vec, error_vec(4,:),'*', test_run_vec, 3*sigma_ydot,'-',test_run_vec,-3*sigma_ydot,'-');
plot(test_run_vec, 3*sigma_ydot,'-');