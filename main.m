% This test script pulls in an example in vivo B1+ field data set along
% with simulated VOPs and a global Q-matrix to allow optimisation for a
% fixed pulse amplitude cardiac imaging sequence to be run.
%
% Created by Arian Beqiri, King's College London, September 2016.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

% Load all data
load bin/b1_data
load bin/Global_Q
load bin/VOPs

N = size(tx,1);     % Size of grid
Nc = size(tx,3);    % Number of coils

nom = 60;           % Nominal flip angle for acquired B1+ maps (degrees)

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints
lSARmax = 10;       % Local SAR constraint (W/kg)
wbSARmax = 4;       % Whole-body SAR constraint (W/kg)

% Sequence Properties
theta_targ = 45;    % Target flip angle (degrees)
delta_1 = 0.538;    % delta_1 as defined in paper
delta_2 = 0.4029;   % delta_2 as defined in paper
t_enc = 1.7;        % Read-out duration (ms)

% Allowed Power Properties
A_amp = 2.5;        % Amplifier constant (W/uT^2)
P_av = 100;         % Average power (W)
P_peak = 1000;      % Peak power (W)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Scale fields and create ROI indices

% Scale transmit field by flip angle
tx = tx/nom;

% Select indices of ROI for heart and over slice
idx = find(M(:));
idx_quad = find(M_quad(:));

%% Target Pattern and A Matrix

% Calculate Quadrature SAR
Quad_Local_SAR = SAR_par(VOP,ones(Nc,1));

% Quad B1 map
quad_map = reshape(sum(tx,3),N*N,1);quad_map(isnan(quad_map))=0;

% Reshape transmit fields
A = reshape(tx, [N*N Nc]);
Amask = A(idx,:);

%% CVX Solution

b1_act_scale = mean(abs(quad_map(idx_quad)))^2;

% Configure variables for optimisation
targmask = ones(length(idx),1);
PA_init = 5;
A_tmpz = sum(A*PA_init,2);
z = exp(1i*(angle(A_tmpz(idx))));
counter = 0;
Bias_conv = [];

% Save variables passed between function calls of optimisation
save('z_tmp','z','counter','Bias_conv')

% Select optimisation to run
opt_sel = questdlg('Select Optimization','Select Optimization',...
                'MSE','Minimum Bias','MSE');

%%% Run Two level Nested Optimisation
switch opt_sel
case 'MSE'
tic
options =optimset('Display','iter','MaxIter',inf,'tolfun',0.1,'tolx',0.01);
[PA_shim,f_val] = fminsearch(@(PA_optim) optim_MSE(PA_optim,A,Amask,VOP,...
QG,targmask,idx,b1_act_scale,Nc,t_enc,lSARmax,wbSARmax,delta_1,delta_2,...
theta_targ,A_amp,P_av,P_peak),PA_init,options);
toc
case 'Minimum Bias'
tic
options =optimset('Display','iter','MaxIter',inf,'tolfun',0.1,'tolx',0.01);
[PA_shim,f_val] = fminsearch(@(PA_optim) optim_minBias(PA_optim,A,Amask,...
VOP,QG,targmask,idx,b1_act_scale,Nc,t_enc,lSARmax,wbSARmax,delta_1,...
delta_2,theta_targ,A_amp,P_av,P_peak),PA_init,options);
toc
end
%%%%%%%%%%%%%%%%%%%%%%

% Optimisation outputs
[rf_dur_shim,t_rms_shim] = SSFP_pulse_params(PA_shim,theta_targ,...
                                                delta_1,delta_2);
TR_shim = max([rf_dur_shim*2 rf_dur_shim+t_enc]);
disp(['Shimmed TR = ' num2str(TR_shim)])

load z_tmp
shim = mls_shim;

%% Calculate Quad with proper power scale and scale to max if max is
% exceeded

% load PAs_TRs_Durs
PAs = 1:0.01:8; TRs = zeros(length(PAs),1);

drmax_f = (sqrt(P_peak/A_amp))./PAs;
dravg_f = zeros(length(PAs),1); T_rms = zeros(length(PAs),1);

parfor ii=1:length(PAs)
    
    [rf_dur_i,T_rms(ii)] = SSFP_pulse_params(PAs(ii),theta_targ,...
                                                delta_1,delta_2);
    TRs(ii) = max([rf_dur_i*2 rf_dur_i+t_enc]);
    
    % Average Power limit
    dravg_f(ii) = sqrt(P_av*TRs(ii)/((PAs(ii).^2*T_rms(ii))*A_amp));
    
end

TR_new = zeros(length(PAs),1); TR_unc = TR_new; Quadrature_Bias = TR_new;
Quadrature_E = TR_new; b1_act_scale_q = TR_new; Quadrature_Var = TR_new;

% Find appropriate PO drive
Quad_PO_drive = 1/mean(abs(quad_map(idx)));

parfor ii=1:length(PAs)

    b1_ampt = PAs(ii);  
    b1_act_scale_q(ii) = (mean(abs(quad_map(idx_quad))))^2*...
                    ((PAs(ii)*Quad_PO_drive)^2)*T_rms(ii)/TRs(ii);

    % Calculate Quad TR to conform to SAR limits
    Quad_SAR = Quad_Local_SAR*b1_act_scale_q(ii);
    TR_new(ii)=round(Quad_SAR/lSARmax*TRs(ii)*100)/100;
    TR_unc(ii)=TR_new(ii);
    if TR_new(ii) < TRs(ii); TR_new(ii) = TRs(ii); end
    
    % Calculate Quadrature Error
    Quadrature_E(ii) = norm((abs(quad_map(idx)*Quad_PO_drive*PAs(ii))-...
        ones(length(idx),1)*b1_ampt))/norm(ones(length(idx),1)*b1_ampt);

    % Calculate Quadrature Bias
    Quadrature_Bias(ii) = abs(b1_ampt-mean(abs(quad_map(idx))*PAs(ii))*...
                            Quad_PO_drive)/b1_ampt*100;
    % Quadrature Variance
    Quadrature_Var(ii) = var(abs(quad_map(idx)*PAs(ii))*Quad_PO_drive,1);
end

figure
plot(TR_unc,PAs,'r','linewidth',2);hold on
plot(TR_new,PAs,'linewidth',2); grid on; axis([0 12 1 12]);
xlabel('TR (ms)'); ylabel('Pulse Amplitude (\muT)')

%% Select Appropriate Quadrature operating point

[TR_quad,ind_quad] = min(TR_new);

b1_act_scale_quad = b1_act_scale_q(ind_quad)*TRs(ind_quad)/TR_quad;
Quadrature_Error = Quadrature_E(ind_quad);

PA_quad = PAs(ind_quad);

%% Final solutions using selected shim
solnFin = squeeze(A)*shim*PA_shim;
Varn = var(abs(solnFin(idx)),1);
solnFin = reshape(solnFin,[N N]);

%% Make Plots and Print SAR Values

Quadrature_SAR = abs([Quad_Local_SAR; ((ones(Nc,1))'*QG*(ones(Nc,1)))])...
   *b1_act_scale_quad;
               
Shimmed_SAR = abs([SAR_par(VOP,shim); (shim'*QG*shim)])*S;

TR_all = [TR_quad TR_shim];
disp(table(Quadrature_SAR,Shimmed_SAR,'RowNames',{'Torso','Whole Body'}));

Shimmed_Error = Error;
Percent_Diff = (Quadrature_Error-Shimmed_Error)/Quadrature_Error*100;

disp(table(Quadrature_Error,Shimmed_Error,Percent_Diff,TR_all))

% Shimmed solution against unshimmed
figure
subplot(121)
imagesc(abs(sum(tx,3))*PA_quad.*M_quad*Quad_PO_drive,...
[0 1.3*PA_quad]);axis off image
title('Quadrature'); h=colorbar('southoutside'); xlabel(h,'\muT')

subplot(122)
imagesc(abs(solnFin).*M_quad,[0 1.3*PA_shim]);axis off image
title('Shimmed'); h=colorbar('southoutside'); xlabel(h,'\muT')

%% Stats
cov_quad = sqrt(Quadrature_Var(ind_quad))/PAs(ind_quad);
cov_shim = sqrt(Varn)/mean(abs(solnFin(idx)));

%% Write out Shims 
fprintf('\nShims = \n\n')
disp(abs(shim))

disp(['Quad TR = ' num2str(TR_quad)])
disp(['Quad Pulse Amplitude = ' num2str(PA_quad)])

disp(['Shim TR = ' num2str(TR_shim)])
disp(['Shim Pulse Amplitude = ' num2str(PA_shim)])