% Mean Squared Error cost function
%
% Created by Arian Beqiri, King's College London, September 2016.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

function [fval] = optim_MSE(P_A,Ap,Amask,VOP,QG,targmask,idx,...
                 b1_act_scale,Nc,t_enc,lSARmax,wbSARmax,delta_1,delta_2,...
                 theta,A_amp,P_av,P_peak)

% Load in MLS phase distribution and counter
zt = load('z_tmp'); z = zt.z; counter = zt.counter;
                
A = Ap * P_A;
AA = Amask * P_A;
T = targmask*P_A;

% Compute pulse duration, RMS duration and TR
[rf_dur,t_rms] = SSFP_pulse_params(P_A,theta,delta_1,delta_2);
TR = max([rf_dur*2 rf_dur+t_enc]);

% SAR Scaling
S = b1_act_scale*(P_A^2)*t_rms/TR;

% Peak and average power drive constraints
drmax = (sqrt(P_peak/A_amp))/P_A;
dravg = sqrt(P_av*TR/((P_A.^2*t_rms)*A_amp));

% Most constraining power constraint
D = min([drmax dravg]);

mlsinit = ones(Nc,1);
mls_shim = zeros(Nc,1);

while (abs(norm(mls_shim - mlsinit)) > 0.1)

    mlsinit = mls_shim;

    %%%% CVX Optimisation %%%%%%%%%%%
    cvx_begin quiet
         variable y(Nc) complex;
         minimise(norm(AA*y - T.*z));
         subject to           
            (y'*QG*y)*S <= wbSARmax;           % Global SAR Constraint
            SAR(y,VOP)*S <= (lSARmax + 0i);    % Local SAR Constraint
            max(abs(y)) <= D;                  % Power constraint
    cvx_end

    mls_shim = y;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mls_tmpsoln = A*mls_shim;  
    z = exp(1i*(angle(mls_tmpsoln(idx))));

    err = norm((abs(mls_tmpsoln(idx))-T))/norm(T);
    conv = abs(norm(mls_shim - mlsinit));
    fprintf('Convergence = %1.6f\tError = %1.4f\n',conv,err);

    % Solutions do not improve much after 1st MLS iteration so force 
    % only one iteration for these
    if counter;break;end
end
disp(['Pulse Amplitude = ' num2str(P_A) '    TR = ' num2str(TR)])

% Compute Bias of solution
Bias = abs(P_A - mean(abs(mls_tmpsoln(idx))))/P_A*100;
Error = err;

% Set counter from zero to one in first function call
counter = 1;

% Plot bias convergence
Bias_conv = [zt.Bias_conv Bias];
plot(Bias_conv);xlabel('Iteration');ylabel('Percentage Bias');grid on
drawnow;hold on
plot([0 length(Bias_conv)+1],[5 5],'--k')

% Cost function value
if Bias <= 5;Bias = 0;end;fval = TR + Bias^2;

save('z_tmp','z','counter','mls_shim','Bias_conv','S','Error')

end