% This test script pulls in an example in vivo B1+ field data set along
% with simulated VOPs and a global Q-matrix to allow optimisation for a
% fixed pulse amplitude cardiac imaging sequence to be run.
%
% Created by Arian Beqiri, King's College London, December 2015.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

% Load all data
load bin/b1_data
load bin/Global_Q
load bin/VOPs
 
N = size(tx,1);     % Size of grid
Nc = size(tx,3);    % Number of coils

%%%%%%%%%%%%%%%%%%%%%% Fixed Pulse and B1 Parameters %%%%%%%%%%%%%%%%%%%%%%

% Nominal flip angle (degrees)
nom = 60;
% B_1^+ pulse amplitude (uT)
b1_amp = 4.5945;
% T_RMS (ms)
t_rms = 0.1764;
% TR (ms)
TR = 6;
% Load Factor from scanner Power Optimisation
LF = 1.78;
% Max drive constraint (20/b1_amp) as defined by RF amplifier limits
drmax = 20/b1_amp;
% Define range of Local SAR constraints in unscaled units
crs = 2.5:-0.4:0.5;
% B1+ based SAR scaling for EM model to quadrature
Bsc = 3.71^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define SAR scaling factor to get SAR in correct units
SAR_scale = Bsc*(b1_amp^2)*t_rms/TR;

% Select indices of ROI for heart and over slice
idx = find(M(:));
idx_quad = find(sum(tx,3));

%% Target Pattern and A Matrix

% Calculate Quadrature values
quad_map_init = reshape(sum(tx,3),N*N,1);
Quad_Local_SAR = SAR_par(VOP,(ones(8,1)));

% Scaling factors based on quadrature normalisation
quad_map = reshape(sum(tx,3),N*N,1);quad_map(isnan(quad_map))=0;

b1_act_scale = mean(abs(quad_map(idx_quad)))^2/(b1_amp^2);
b1_act_scale_quad = b1_act_scale*(LF^2);

% Define target pattern and matrix of coil sensitivities
targ = ones(N*N,1)*b1_amp;
A = reshape(tx, [N*N Nc]);

% Set indices of ROI and set A and target against these
targmask = targ(idx);
Amask = A(idx,:);

% Calculate Quadrature Error in ROI
Quadrature_Error = norm((abs(quad_map(idx)*LF)-targ(idx)))/norm(targ(idx));

% Calculate Quad with proper power scale and scale to max if max is
% exceeded
Quad_PO_drive = (0.95*b1_amp)/mean(abs(quad_map(idx)));
if Quad_PO_drive > drmax; Quad_PO_drive = drmax;end

Quadrature_Error_PO = norm((abs(quad_map(idx)*Quad_PO_drive)-targ(idx)))/norm(targ(idx));
b1_act_scale_quad_PO = mean(abs(quad_map(idx_quad)*Quad_PO_drive))^2/(b1_amp^2);

%% CVX Solution

% Calculate CVX L-curve
[all_shims,Error,LSAR,GSAR] = CVX_b1_shim_LCurve(A,Amask,VOP,QG,targmask,idx,drmax,crs,Nc);

% Plot L-curve
figure
subplot(121)
scatter(Quadrature_Error,Quad_Local_SAR*SAR_scale*b1_act_scale_quad)
hold on
scatter(Quadrature_Error_PO,Quad_Local_SAR*SAR_scale*b1_act_scale_quad_PO)
hold on
plot(Error,LSAR*SAR_scale*b1_act_scale)
title('Max Local SAR')
subplot(122)
scatter(Quadrature_Error,abs(((ones(8,1)))'*QG*((ones(8,1))))*SAR_scale*b1_act_scale_quad)
hold on
scatter(Quadrature_Error_PO,abs(((ones(8,1)))'*QG*((ones(8,1))))*SAR_scale*b1_act_scale_quad_PO)
hold on
plot(Error,GSAR*SAR_scale*b1_act_scale)
title('Global SAR')

% Select shim that achieves at least 95% of FA
tmpsoln = Amask*(all_shims.');
ind_lambda = find(mean(abs(tmpsoln))/b1_amp>0.95);
shim = (all_shims(ind_lambda(end),:)).';

% Final solutions using selected shim
solnFin = A*shim;
solnFin = reshape(solnFin,[N N]);

%% Make Plots and Print SAR Values
Quadrature_SAR = abs([Quad_Local_SAR;...
            ((ones(8,1))'*QG*(ones(8,1)))])*SAR_scale*b1_act_scale_quad;

Quadrature_SAR_PO = abs([Quad_Local_SAR;...
            ((ones(8,1))'*QG*(ones(8,1)))])*SAR_scale*b1_act_scale_quad_PO;
               
Shimmed_SAR = abs([SAR_par(VOP,(shim));...
                ((shim)'*QG*(shim))])*SAR_scale*b1_act_scale;

table(Quadrature_SAR,Quadrature_SAR_PO,Shimmed_SAR,'RowNames',{'Local','Global'})

Shimmed_Error = Error(ind_lambda(end));
Percent_Diff = (Quadrature_Error-Shimmed_Error)/Quadrature_Error*100;
Percent_Diff_PO = (Quadrature_Error_PO-Shimmed_Error)/Quadrature_Error_PO*100;

table(Quadrature_Error,Quadrature_Error_PO,Shimmed_Error,Percent_Diff,Percent_Diff_PO)

SAR_Percent_Diff_PO = (Quadrature_SAR_PO(1) - Shimmed_SAR(1))/Quadrature_SAR_PO(1)*100;

table(Percent_Diff_PO,SAR_Percent_Diff_PO)

% Shimmed solution against unshimmed
figure
scale = [0 1.4*b1_amp];
subplot(121)
imagesc(abs(sum(tx,3)*Quad_PO_drive),scale);axis off
title('Quadrature')
colorbar('southoutside')

subplot(122)
imagesc(abs(solnFin),scale);axis off
title('CVX Shim')
colorbar('southoutside')

% Display Shims 
fprintf('\nShim amplitudes = \n\n')
disp(abs(shim))