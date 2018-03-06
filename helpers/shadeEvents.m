function h = shadeEvents(ax,t,ev,varargin)
% shadeEvents(ax,t,ev,varargin)
%
% Shade regions on an axis, 'ax', when an event is occurring. To do this
% supply the logical vector of the events and the corresponding time
% vector. Plotting parameters can be supplied as well.

if ~exist('varargin','var') || length(varargin) < 1
    varargin = {[0.8 0.8 0.8], 'FaceAlpha',0.6,'linestyle','none'};
end

evOn = diff(ev) > 0; % event onset
evOff = diff(ev) < 0; % event offset
nEvs = sum(evOn); % how many events

if ~isrow(t); t = t'; end % make it a column vector
tEvs = reshape([t(evOn); t(evOff)],1,[]);

if nargout > 1
    h = fill(ax,repelem(tEvs,2),repmat([ax.YLim, fliplr(ax.YLim)], nEvs),varargin{:});
else
    fill(ax,repelem(tEvs,2),repmat([ax.YLim, fliplr(ax.YLim)], nEvs),varargin{:});
end

end