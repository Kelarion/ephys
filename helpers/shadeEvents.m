function h = shadeEvents(ev,col,t,varargin)
% shadeEvents(ev,t,ax[,Name,Value])
%
% Shade regions on an axis, 'ax', when an event is occurring. 
%
% Inputs: 
%   - ev: either [Onset Offset] matrix or logical vector
%   - col: color, passed straight into 'fill' function
%   - t: if ev is logical, supply the corresponding time vector
%        (otherwise supply empty array)
% Extra arguments are passed into 'fill' function.

ax = gca;
if ~exist('col','var') || isempty(col), col=[0.8 0.8 0.8]; end

if ~exist('varargin','var') || length(varargin) < 1
    varargin = {'FaceAlpha',0.3,'linestyle','none'};
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
    h = fill(ax,repelem(tEvs,2),repmat([ax.YLim, fliplr(ax.YLim)], nEvs),col,varargin{:});
else
    fill(ax,repelem(tEvs,2),repmat([ax.YLim, fliplr(ax.YLim)], nEvs),col,varargin{:});
end

end