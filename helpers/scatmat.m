function Ax = scatmat(X,ptype,varargin)
% Ax = scatmat(X,ptype[,Name,Value])
%
% X is [nObs,nVar]
% ptype is either 'scatter' (default) or 'density'
% specify axes labels as (Name,Value) pair
% Ax is [nVar,nVar]

if nargin<2, ptype = 'scatter'; end

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
switch ptype % a bit more flexible than just using 'plotmatrix' function
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
                if row == col
                    Ax(row,col) = subplot(nVar,nVar,j);
                    histogram(X(:,row),nbin,'edgealpha',0);
                    theseXlims = xlim;
                else
                    Ax(row,col) = subplot(nVar,nVar,j);
                    
                    [n, c] = hist3(X(:,[col, row]),[nbin nbin]); % make selection of nbin smarter
                    imagesc(c{1},flipud(c{2}(:)),flipud(n')); colormap gray;
                    set(Ax(row,col),'Ydir','normal')
                end
                j = j+1;
            end
        end
    otherwise
        error('Please supply a valid plot type, or don''t supply one at all')
end
if ~isempty(labs) % label axes
    for v = 1:nVar
        ylabel(Ax(v,1),labs{v})
        xlabel(Ax(nVar,v),labs{v})
    end
end

end