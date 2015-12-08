function [SAR] = SAR_par(VOPm,x)

SAR = zeros(size(VOPm,1),1);

parfor k=1:size(VOPm,1)
    SAR(k) = (x'*(squeeze(VOPm(k,:,:)))*x);
end

[SAR] = max(SAR);