function [val,mag] = circmean(X,dim,nanflag)
% [val,mag] = circmean(X[,dim,nanflag])
%
% Calculate the mean of points that lie on a circle
% Takes mean along rows by default

if isvector(X), X = X(:); end % make column

if nargin<2, dim = 1; end % set defaults
if nargin<3, nanflag = false; end

% place along unit circle, take mean, get phase
if nanflag 
    val = angle(nanmean(exp(1i*X),dim));
    mag = abs(nanmean(exp(1i*X),dim));
else
    val = angle(mean(exp(1i*X),dim));
    mag = abs(mean(exp(1i*X),dim));
end