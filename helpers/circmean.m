function [val,mag] = circmean(X,dim,nanflag)
% [val,mag] = circmean(X[,dim,nanflag])
%
% Calculate the mean of points that lie on a circle
% Takes mean along longest dimension by default

if nargin<2, [~,dim] = max(size(X)); end
if nargin<3, nanflag = false; end

% place along unit circle, take mean, get phase
if nanflag 
    val = angle(nanmean(exp(1i*X),dim));
    mag = abs(nanmean(exp(1i*X),dim));
else
    val = angle(mean(exp(1i*X),dim));
    mag = abs(mean(exp(1i*X),dim));
end