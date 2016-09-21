% Function to output RF pulse timings given flip angle, pulse amplitude,
% delta_1 and delta_2
%
% Created by Arian Beqiri, King's College London, September 2016.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

function [rf_dur,t_rms] = SSFP_pulse_params(b1_amp,theta,delta_1,delta_2)

% Convert flip angle to radians in radians
theta_rad = theta/180*pi;

% Define gamma in unit of rad s^-3 uT^-1
gamma = 0.2675;

% RF pulse duration (ms)
rf_dur = theta_rad/(gamma*b1_amp*delta_1);

% RF pulse RMS time (ms)
t_rms = delta_2 * rf_dur;

end