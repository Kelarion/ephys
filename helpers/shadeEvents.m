function h = shadeEvents(ax,ev,t,varargin)
% shadeEvents(ax,ev,t[,Name,Value])
%
% Shade regions on an axis, 'ax', when an event is occurring. 
%
% Inputs: 
%   - ax: axes
%   - ev: either [Onset Offset] matrix or logical vector
%   - t: if ev is logical, supply the corresponding time vector
%        (otherwise supply empty array)
% Plotting parameters can be supplied as well.

if ~exist('varargin','var') || length(varargin) < 1
    varargin = {[0.8 0.8 0.8], 'FaceAlpha',0.6,'linestyle','none'};
end

if (any(ev) > 1)
    t = t(:);
    evOn = diff([0 ev]) > 0; % event onset
    evOff = diff([ev 0]) < 0; % event offset
    tEvs = reshape([t(evOn); t(evOff)],1,[]);
else
    if diff(size(ev))>0
        ev = ev';
    end
    evOn = ev(:,1);
    evOff = ev(:,2);
    tEvs = reshape([evOn; evOff],1,[]);
end
nEvs = length(tEvs)/2; % how many events

if nargout > 1
    h = fill(ax,repelem(tEvs,2),repmat([ax.YLim, fliplr(ax.YLim)], nEvs),varargin{:});
else
    fill(ax,repelem(tEvs,2),repmat([ax.YLim, fliplr(ax.YLim)], nEvs),varargin{:});
end

end