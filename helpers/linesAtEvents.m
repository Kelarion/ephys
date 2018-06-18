function linesAtEvents(tEvs,ax,varargin)
% linesAtEvents(tEvs,ax[,Name,Value])
% 
% Draws a vertical line at all events in tEvs. Remaining args are passed
% into plot directly. Make sure the axes ylimits are set.

if ~exist('ax','var') || isempty(ax)
    ax = gca;
end

nEvs = length(tEvs);
plot(ax,repelem(tEvs,3),repmat([ax.YLim NaN], nEvs),varargin{:})