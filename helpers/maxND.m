function [val, inds] = maxND(A)
% [val inds] = maxNd(A)

dmnlty = sum(size(A) > 1);

inds = zeros(dmnlty,1);
vals = zeros(dmnlty,1);
for iDim = 1:dmnlty
    [vali, indi] = max(max(A,[],dmnlty - iDim + 1));
    inds(iDim) = indi;
    vals(iDim) = vali;
end
val = mean(vals);

end