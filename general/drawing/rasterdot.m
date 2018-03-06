function order = rasterdot(st, varargin)
%   Code originally by original by Ming Gu, changes by CMW and AWH.
%
%RASTER   Spike raster plot.
%
%   RASTER(ST) given a cell array of spike times ST generates a raster plot.
%
%   Optional args:
%   'axes'        Axes handle to plot into - gca if empty
%   'color'       Colormap matrix
%   'doblack'     Set figure to "black" colordef
%   'line_width'  Set linewidth value
%   'spacing'     Spacing of tick marks on y axis (default: 1)
%   'tick_height' Height of the tick mark (default: 0.75)
%   'order'       Row order in which cells should be displayed
%   'tlim'        Time limits or empty for full range
%   'byrate'      Order cells by number of spikes; this will
%                 override 'order'. (default: false)
%   'use_3d'      Use 3D plot instead (default: false)
%

if isnumeric(st)
    st = {st};
end

ip = inputParser();
ip.addParamValue('axes', [], @(x) isempty(x) || ishandle(x));
ip.addParamValue('color', lines(7), @(x) isnumeric(x) && (size(x,2)==3));
ip.addParamValue('doblack', false, ...
    @(x) isscalar(x) && (islogical(x) || isnumeric(x)));
ip.addParamValue('line_width', 1, @(x) isscalar(x) && isnumeric(x));
ip.addParamValue('spacing', 1, @(x) isscalar(x) && isnumeric(x));
ip.addParamValue('tick_height', .75, @(x) isscalar(x) && isnumeric(x));
ip.addParamValue('order', 1:numel(st), ...
    @(x) isnumeric(x) && (numel(x) == numel(st)));
ip.addParamValue('tlim', [], ...
    @(x) isempty(x) || (isnumeric(x) && (numel(x) == 2)));
ip.addParamValue('byrate', false, ...
    @(x) isscalar(x) && (isnumeric(x) || islogical(x)));
ip.addParamValue('use_3d', false, ...
    @(x) isscalar(x) && (isnumeric(x) || islogical(x)));
ip.parse(varargin{:});
params = ip.Results;


if isempty(params.axes)
    params.axes = gca;
end

if params.doblack
    fh = ancestor(params.axes, 'figure');
    colordef(fh, 'black');
end

if params.byrate
    [~, params.order] = sort( cellfun(@numel, st), 'descend' );
end


nclr = size(params.color, 1);
st = st(params.order);
ax = params.axes;

% restrict to time limits
if ~isempty(params.tlim)
    f = @(x,y)(x(x >= params.tlim(1) & x <= params.tlim(2)));
    st = cellfun(f, st, 'UniformOutput', false);
end
sel_nonempty = find(~cellfun(@isempty, st));

for ii = sel_nonempty(:)'
    clrid = mod(params.order(ii) - 1, nclr) + 1;
    clrid = mod(ii - 1, nclr) + 1;
    nt = length(st{ii});
    u = repmat(st{ii}(:)', 3, 1);
    v = repmat([1-params.tick_height/2; 1+params.tick_height/2; nan], 1, nt) + ...
        params.spacing*(ii-1);
    if params.use_3d
        plot3(ax, u(:), v(:), v(:)*0,  '-', 'color', params.color(clrid,:), ...
            'LineWidth', params.line_width);
        view(ax, 2);
    else
        plot(ax, u(:), v(:), '-', 'color', params.color(clrid,:), ...
            'LineWidth', params.line_width);
    end
    hold(ax, 'on')
end
% axis(ax, 'ij')
if ~isempty(params.tlim)
    xlim(ax, params.tlim);
end

if ~isempty( ii )
    ylim(ax, params.spacing * [0 numel(st)] + 0.5);
end
hold(ax, 'off')

if params.doblack
    colordef(fh, 'white');
end

order = params.order;
