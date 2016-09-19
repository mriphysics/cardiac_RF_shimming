% Function to find max local SAR for a given shim.
%
% Created by Arian Beqiri, King's College London, December 2015.
% Email: arian.beqiri@kcl.ac.uk
%
% This code is free under the terms of the MIT license.

function [SAR] = SAR(x,Q)

% Compute all local SAR values
for k=1:size(Q,1)
    SAR(k) = ((x)'*(squeeze(Q(k,:,:)))*(x));
end

% Output max local SAR
[SAR] = max(SAR);
end
