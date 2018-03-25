function [inds, nearest] = findNearestPoint(toThese, inThese)
% function [inds, nearest] = findNearestPoint(toThese, inThese)
% 
% Finds the value of "inThese" closest to a value of "toThese", and the
% index of it, for all "toThese".
%
% inThese and ofThese must both already be sorted. 
%
% modified from Nick Steinmetz (I think), CortexLab, UCL 

toThese = toThese(:)';% to row
inThese = inThese(:)'; 

[~,ii] = sort([inThese toThese]);
[~,ii2] = sort(ii); 

prevInds = ii2(length(inThese)+1:end)-1-(0:length(toThese)-1);

nearest = zeros(size(toThese));
inds = zeros(size(toThese));

tooEarly = prevInds<1;
nearest(tooEarly) = inThese(1);
inds(tooEarly) = 1;

tooLate = prevInds>=length(inThese);
nearest(tooLate) = inThese(end);
inds(tooLate) = length(inThese);

nextDiff = abs(toThese(~tooEarly&~tooLate)-inThese(prevInds(~tooEarly&~tooLate)+1));
prevDiff = abs(toThese(~tooEarly&~tooLate)-inThese(prevInds(~tooEarly&~tooLate)));

nextClosest = false(size(toThese));
prevClosest = false(size(toThese));
nextClosest(~tooEarly&~tooLate) = nextDiff<prevDiff;
nearest(nextClosest) = inThese(nextClosest);
inds(nextClosest) = prevInds(nextClosest)+1;

prevClosest(~tooEarly&~tooLate) = nextDiff>=prevDiff;
nearest(prevClosest) = inThese(prevClosest);
inds(prevClosest) = prevInds(prevClosest);