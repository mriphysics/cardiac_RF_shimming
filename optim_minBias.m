% Minimum Bias cost function
%
% Created by Arian Beqiri, King's College London, September 2016.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

function [fval] = optim_minBias(P_A,Ap,Amask,...
                        VOP,QG,targmask,idx,b1_act_scale,Nc,t_enc)

% Load in MLS phase distribution and counter
zt = load('z_tmp'); z = zt.z; counter = zt.counter;
                
A = Ap * P_A;
AA = Amask * P_A;
T = P_A;

% Compute pulse duration, RMS duration and TR
[~,rf_dur,t_rms] = SSFP_calc_duty_cycle(P_A);
TR = max([rf_dur*2 rf_dur+t_enc]);

% SAR Scaling
S = b1_act_scale*(P_A^2)*t_rms/TR;

% Peak and average power drive constraints
drmax = 20/P_A;
dravg = sqrt(100*TR/((P_A.^2*t_rms)*2.25));

% Most constraining power constraint
D = min([drmax dravg]);

sinit = ones(Nc,1);
mls_shim = zeros(Nc,1);

len_z = length(z);

while (abs(norm(mls_shim - sinit)) > 0.1)

    sinit = mls_shim;

    %%% CVX Optimisation %%%%%%%%%%%
    cvx_begin quiet
         variable y(Nc) complex;
         minimise(norm(sum((AA*y).*(conj(z)))/len_z - T));
         subject to           
            (y'*QG*y)*S <= 4;           % Global SAR Constraint
            SAR(y,VOP)*S <= (10 + 0i);  % Local SAR Constraint
            max(abs(y)) <= D;           % Power Constraint
    cvx_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mls_shim = y;
    mls_tmpsoln = A*mls_shim;  
    z = exp(1i*(angle(mls_tmpsoln(idx))));

    err = norm((abs(mls_tmpsoln(idx))-T*targmask))/norm(T*targmask);
    conv = abs(norm(mls_shim - sinit));
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
plot(Bias_conv)
drawnow;hold on

% Cost function value
fval = TR + Bias^2;

save('z_tmp','z','counter','mls_shim','Bias_conv','S','Error')

end