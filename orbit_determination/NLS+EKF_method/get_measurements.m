function [y] = get_measurements(time, fem_states, mom_states, freq, radius_planet, flag_for_test)
%GETMEASUREMENTS Generate the synthetic Doppler shift measurements
%   Given the states of the mothersat and femtosat, generate
%   synthetic Doppler shift measurements.
% 
%   Inputs:
%     time          : array of time instances for fem and mom states
%     fem_states    : matrix of the femtosatellite states
%     mom_states    : matrix of the mothersatellite's states
%     freq          : transmission frequency
%     radius_planet : radius of the planetary body being orbitted
%     flag_for_test : 1 = test to see if line-of-sight is obstructed
%                     0 = skip test // see footnote [1]
%
%   Output:
%     y             : 2-column matrix, [time measurements]
%

m             = size(fem_states, 1);         % number of measurements
n             = size(fem_states, 2);         % number of states
c             = 3e8;                         % speed of light, m/s
lamda         = (c/freq)/1000;               % wavelength, km
y             = zeros(0,2);                  % create empty array
meas_index    = 1;

for ii = 1:m
  x_fem       =  fem_states(ii,1);
  y_fem       =  fem_states(ii,2);
  x_mom       =  mom_states(ii,1);
  y_mom       =  mom_states(ii,2);
  r           =  [x_fem; y_fem];
  R           =  [x_mom; y_mom];
  r_dot       =  [fem_states(ii,3); fem_states(ii,4)];
  R_dot       =  [mom_states(ii,3); mom_states(ii,4)];
  rho_dot_mag =  get_rho_dot_mag(r,r_dot,R,R_dot);
  angle_mom   =  atan2(y_mom, x_mom);
  angle_fem   =  atan2(y_fem, x_fem);
  temp        =  asin(radius_planet/R);
  angle_test  =  pi/2 - temp;
  max_angle   =  angle_mom + angle_test;
  min_angle   =  angle_mom - angle_test;

  % Here, we need to test for line-of-sight. There's a two-part test.
  % Also, we need to take care of the case when the two sats have the
  % same y-value in the coordinate system.
  slope  =  (y_fem-y_mom)/(x_fem-x_mom);
  if slope == 0
    % The sats have the same y-value on the coordinate frame.
    if flag_for_test
      % Use the angle test
      if (angle_fem < max_angle && angle_fem > min_angle)
        % If they can see each other, write out measurement
        y(meas_index, :) = [time(ii) -rho_dot_mag/lamda];
        meas_index       = meas_index + 1;
      end % angle test
    else
      % no test for obstruction, write out measurement
      y(meas_index, :) = [time(ii) -rho_dot_mag/lamda]; % Dopp. shift
      meas_index       = meas_index + 1;      
    end % if (flag_for_test) && (abs(y_fem) > radius_planet)
  else
    b            =  y_fem - slope*x_fem;
    slope_perpen = -1/slope;
    x_intersect  =  b/(slope_perpen - slope);
    y_intersect  =  slope*x_intersect + b;
    rad_perp     =  norm([x_intersect y_intersect]);

    if flag_for_test
      if (angle_fem < max_angle && angle_fem > min_angle) || ...
         (rad_perp > radius_planet)
        % line-of-sight not obstructed by planetary body   
        y(meas_index, :) = [time(ii) -rho_dot_mag/lamda];
        meas_index       = meas_index + 1;
      end
    else
      % write all measurements out
      y(meas_index, :) = [time(ii) -rho_dot_mag/lamda]; % Dopp. shift
      meas_index       = meas_index + 1;
    end % if flag_for_test
  end % if slope == 0
end % for ii = 1:m

% [1] Reason for conditional test:
%   When creating measurements for the simulation, we want to test if
% the line-of-sight is obstructed by the planet, to understand when
% measurements would and wouldn't be available.
%   However, in the Kalman filter loop, we're using the state
% estimate. So, during the loop, we may have a measurement for an
% instant, but because of the femtosat state estimate & geometry,
% this function would say that no measurement is available, and we'd
% have to ignore the measurement data point during that iteration
% of the loop. Rather than ignore measurement data, we'll just ignore
% the geometry of the planet blocking the signal.