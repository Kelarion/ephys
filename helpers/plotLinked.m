function viewTogether(sig1, sig2, blackBackground)
% viewTogether(sig1, sig2, blackBackground)
%
% Don't use this function programmatically.

if nargin < 3
    blackBackground = 1;
end
fig = figure('Units','Normalized','Position', [0.1495 0.2306 0.7828 0.6907]); 
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);
plot(ax1,  sig1)
plot(ax2,  sig2)

linkaxes([ax1, ax2], 'x');
if blackBackground
    whitebg(fig,'k');
end

end