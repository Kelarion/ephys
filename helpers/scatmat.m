function Ax = scatmat(X,ptype,varargin)
% Ax = scatmat(X,ptype[,Name,Value])
%
% Scatter-plot matrix with histograms along the diagonal, with the option
% to plot density rather than points and specify axis labels.
% 
% Inputs: 
%  - X is [nObs,nVar]
% Optional:
%  - ptype is either 'scatter' (default) or 'density'
%  - Specify axes labels as (X,__,'labels', {lab1,lab2,...})
%  - (for density plot) Specify the number of bins as (X,__,'nbin',val)
%  - (for scatter plot) Aditional arguments for the 'plotmatrix' function
% Outputs:
%  - Ax is [nVar,nVar]

if ~exist('ptype','var'), ptype = 'scatter'; end

if ~isempty(varargin) % parse varargin
    nms = lower(varargin(1:2:end));
    rgs = varargin(2:2:end);
    include = false(1,length(varargin)/2);
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
    nExtraArgs = length(varargin) - 2*sum(include);
    extraArgs = cell(nExtraArgs,1);
    extraArgs(1:2:end) = nms(~include);
    extraArgs(2:2:end) = rgs(~include);
else
    nbin = 15;
    labs = {};
    extraArgs = {};
end

nVar = size(X,2);
switch ptype % a bit more flexibility than the 'plotmatrix' function
    case 'scatter'
        if isempty(extraArgs)
            [~,Ax,~,H] = plotmatrix(X);
        else
            [~,Ax,~,H] = plotmatrix(X,extraArgs{:});
        end
        for jj = 1:length(H) % make histograms prettier
            H(jj).EdgeAlpha = 0;
            H(jj).FaceColor = 'r';
        end
    case 'density'
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
                    imagesc(c{1},flipud(c{2}(:)),flipud(n')); colormap gray;
                    set(Ax(row,col),'Ydir','normal')
                    
                    if row<nVar, Ax(row,col).XTick = []; end
                    if col>1, Ax(row,col).YTick = []; end
                end
                j = j+1;
            end
        end
    otherwise
        error('Please supply a valid plot type, or don''t supply one at all')
end
if ~isempty(labs) % label axes
    for v = 1:nVar
        ylabel(Ax(v,1),labs{v},'Rotation',60,'horizontalalignment','right')
        xlabel(Ax(nVar,v),labs{v})
    end
end

end