function [duty,rf_dur,t_rms] = SSFP_calc_duty_cycle(b1_amp)

%%%%%%%%%% TFE Default Sequence Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RF sequence FA in radians
theta = 45/180*pi;
% B_1^+ pulse amplitude (uT)
b1_amp_default = 3.6915;
% T_RMS (ms)
t_rms_default = 0.5957;
% RF duration (ms)
rf_dur_default = 1.4784;
% TR (ms)
TR = 3.24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate ct_rms from defaults
ct_rms = t_rms_default/rf_dur_default;

% Define gamma in unit of rad s^-3 uT^-1
gamma = 0.2675;

% Calculate t_eff
ct_eff = theta/(gamma*b1_amp_default*rf_dur_default);

% Calcute duty cycle
duty = (theta/b1_amp)*(ct_rms/(gamma*ct_eff*TR));

rf_dur = theta/(gamma*b1_amp*ct_eff);

t_rms = ct_rms * rf_dur;

end
