function Ax = scatmat(X,ptype,varargin)
% Ax = scatmat(X,ptype[,Name,Value])
%
% Scatter-plot matrix with histograms along the diagonal, with the option
% to plot density rather than points and to specify axis labels.
% 
% Inputs: 
%  - X is [nObs,nVar]
% Optional:
%  - ptype is either 'scatter' (default) or 'density'
%  - Specify axes labels as (X,__,'labels', {lab1,lab2,...})
%  - Specify the number of histogram/density bins as (X,__,'nbin',val)
%  - (for scatter plot) Aditional arguments for the 'scatter' function
% Output (optional):
%  - Ax is [nVar,nVar] axes

if ~exist('ptype','var'), ptype = 'scatter'; end

[nbin,labs,marker,extraArgs] = parseArgs(varargin);

gcf;

nVar = size(X,2);
switch ptype 
    case 'scatter' % more flexibility than the 'plotmatrix' function
        j = 1;
        for row = 1:nVar
            for col = 1:nVar
                Ax(row,col) = subtightplot(nVar,nVar,j,[0.005 0.005]);
                varXlim = [min(X(:,col)) max(X(:,col))]; % set lims for all plots
                varYlim = [min(X(:,row)) max(X(:,row))];
                rowCntr = linspace(varXlim(1),varXlim(2),nbin);
                if row == col
                    d1 = mean(diff(rowCntr))/2;
                    histogram(X(:,row),rowCntr-d1,'facecolor','b','edgealpha',0);
                    
                    if row<nVar, Ax(row,col).XTick = []; end
                    if col>1, Ax(row,col).YTick = []; end
                    Ax(row,col).XLim = varXlim;
                else
                    scatter(X(:,col),X(:,row),marker,extraArgs{:});
                    
                    if row<nVar, Ax(row,col).XTick = []; end
                    if col>1, Ax(row,col).YTick = []; end
                    Ax(row,col).XLim = varXlim;
                    Ax(row,col).YLim = varYlim;
                end
                j = j+1;
            end
        end
    case 'density' % needs some work
        j = 1;
        for row = 1:nVar
            for col = 1:nVar
                Ax(row,col) = subtightplot(nVar,nVar,j,[0.005 0.005]);
                varXlim = [min(X(:,col)) max(X(:,col))];
                varYlim = [min(X(:,row)) max(X(:,row))];
                rowCntr = linspace(varXlim(1),varXlim(2),nbin);
                colCntr = linspace(varYlim(1),varYlim(2),nbin);
                if row == col
                    d1 = mean(diff(rowCntr))/2;
                    histogram(X(:,row),rowCntr-d1,'edgealpha',0);
                    
                    if row<nVar, Ax(row,col).XTick = []; end
                    if col>1, Ax(row,col).YTick = []; end
                    Ax(row,col).XLim = varXlim;
                else
                    [n, c] = hist3(X(:,[col, row]),{rowCntr colCntr});
                    imagesc(c{1},flipud(c{2}(:)),flipud(n')); 
                    colormap gray;
                    set(Ax(row,col),'Ydir','normal')
                    q = abs(quantile(n(:),0.85));
                    caxis([0 q])
                    
                    if row<nVar, Ax(row,col).XTick = []; end
                    if col>1, Ax(row,col).YTick = []; end
                end
                j = j+1;
            end
        end
    otherwise
        error(['''' ptype ''' not a valid plot type. Please supply a valid plot type, or don''t supply one at all'])
end
if ~isempty(labs) % label axes
    for v = 1:nVar
        ylabel(Ax(v,1),labs{v},'Rotation',60,'horizontalalignment','right')
        xlabel(Ax(nVar,v),labs{v})
    end
end
for v = 1:nVar
    linkaxes(Ax(:,v),'x')
    linkaxes(Ax(v,setdiff(1:nVar,v)),'y')
end
linkaxes(Ax(1,2:end),'y')

end

%% arg parser
function [nbin,labs,marker,extraArgs] = parseArgs(args)

if ~isempty(args) % parse varargin
    nms = lower(args(1:2:end));
    rgs = args(2:2:end);
    include = false(1,length(args)/2);
    if any(contains(nms,'nbin'))
        include = include | contains(nms,'nbin');
        nbin = rgs{contains(nms,'nbin')};
    else
        nbin = 15;
    end
    if any(contains(nms,'labels'))
        include = include | contains(nms,'labels');
        labs = rgs{contains(nms,'labels')};
    else
        labs = {};
    end
    if any(contains(nms,'mkr'))
        include = include | contains(nms,'mkr');
        marker = rgs{contains(nms,'mkr')};
    else
        marker = '.';
    end
    nExtraArgs = length(args) - 2*sum(include);
    extraArgs = cell(nExtraArgs,1);
    extraArgs(1:2:end) = nms(~include);
    extraArgs(2:2:end) = rgs(~include);
else
    nbin = 15;
    labs = {};
    marker = '.';
    extraArgs = {};
end

end

%% helper
function lims = limsAround(cntr,data) % tried this for more robust limits, but doesn't really help
% lims = limsAround(cntr,data)
%
% make limits aroud data that are centered around cntr

lims = [min(data) max(data)];
[~,movelim] = max(abs(cntr-fliplr(lims))); % flipped for ease of indexing
d = 2*cntr - (sum(lims));
lims(movelim) =  lims(movelim) + d;


end
