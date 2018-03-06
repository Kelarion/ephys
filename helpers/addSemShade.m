function h = addSemShade(ax,ex,why,dim,varargin)
% h = addSemShade(ax,ex,why,dim,varargin)
%
% specify the dimension to compute stats over as 'dim'. varargin are the
% plotting parameters (supplied to 'fill').

if ~exist('dim','var') || isempty(dim)
    dim = 1;
end
if ~exist('varargin','var') || length(varargin) < 1
    varargin = {[0.8 0.8 0.8], 'FaceAlpha',0.7,'linestyle','none'};
end


meen = nanmean(why,dim); seem = nanstd(why,[],dim)/sqrt(size(why,dim));

% if size(why,dim) <= 2
%     seem = rand(size(meen)).*meen*2;
% end

h = fill(ax,[ex,flipud(ex')'],[(meen+seem),flipud([meen-seem]')'],varargin{:});

end

