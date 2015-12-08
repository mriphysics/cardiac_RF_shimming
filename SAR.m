%% Function to find max local SAR for a given shim

function [SAR] = SAR(x,Q)

% Compute all local SAR values
for k=1:size(Q,1)
    SAR(k) = ((x)'*(squeeze(Q(k,:,:)))*(x));
end

% Output max local SAR
[SAR] = max(SAR);
end
