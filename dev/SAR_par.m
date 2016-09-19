% Parallelised function to find max local SAR for a given shim which can be
% used outside of CVX.
%
% Created by Arian Beqiri, King's College London, December 2015.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

function [SAR] = SAR_par(VOPm,x)

SAR = zeros(size(VOPm,1),1);

parfor k=1:size(VOPm,1)
    SAR(k) = (x'*(squeeze(VOPm(k,:,:)))*x);
end

[SAR] = max(SAR);