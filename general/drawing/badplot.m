function badplot(theseData,X,t,labs,useCursor,varargin)
% badplot(theseData,X[,t,labels,useCursor])
% 
% Plot a function (X) of the datapoint at the current mouse position. Works
% for both continuous position and data cursor position.
%
% ATM, X needs to be a matrix where column corresponds to a dataponit in
% theseData, or a cell vector where each entry corresponds to a datapoint.

supportedTypes = {'Scatter','Line'};
if nargin<4, useCursor = false; end
if nargin<5, labs = []; end
labs = split(num2str(labs(:)'));

% get initial settings
if iscell(X)
    nPnts = length(X);
else
    nPnts = size(X,2);
end

switch whichClass(theseData)
    case 'Axes' % more complicated, they only gave the axes
        targetAx = theseData;
        % see which data in target axis might be the right one
        dats = allchild(targetAx);
        scats = dats(contains(whichClass(dats),supportedTypes));
        if isempty(scats)
            msg = ['No plots in given axes are of supported type ('  ...
                strjoin(supportedTypes,', ') ')'];
            error(msg)
        else
            candidateData = [];
            for ii = 1:length(scats)
                lx = length(scats(ii).XData);
                ly = length(scats(ii).YData);
                lz = length(scats(ii).ZData);
                if (lx == nPnts) && (ly == lx) && ((ly == lz) || ~lz)
                    candidateData = [candidateData scats(ii)];
                end
            end
            if isempty(candidateData)
                msg = ['No plots in given axes of supported type ('  ...
                    strjoin(supportedTypes,', ') ') match data in ''X'''];
                error(msg)
            elseif length(candidateData)>1
                warning('Multiple data in axes match input, chosing one')
                whichDat = randi(length(candidateData),1);
                targetData = candidateData(whichDat);
            else 
                targetData = candidateData;
            end
        end
        
        if eq(targetAx,gca)
            figure;
            ax = axes();
        else
            ax = gca;
        end
    case supportedTypes % easier, they gave the plot handle
        targetData = theseData;
        targetAx = get(theseData,'parent');
        if eq(get(theseData,'parent'),gca)
            figure;
            ax = axes();
        else
            ax = gca;
        end
    otherwise
        msg = ['Data specified are not of a supported a class ('  ...
            strjoin(supportedTypes,', ') ')'];
        error(msg)
end

targetCoords = [targetData.XData;targetData.YData;targetData.ZData]';
if (nargin<3) && (~iscell(X))
    t = 1:size(X,1); % for now
end

if iscell(X)
    dataPlot = imagesc(ax,nan(size(X{1})));
else
    dataPlot = plot(ax,t,nan(1,size(X,1)));
end

fig = get(targetAx,'parent');
tabulaRasa(fig);
if useCursor
    set(fig,'units','normalized','WindowButtonMotionFcn',{@takeBreath,targetAx.Position})
    dcm_obj = datacursormode(fig);
    lh = addlistener(fig,'CurrentPoint','PostSet', ...
        @(src,evnt)updateTextCsr(src,evnt,dataPlot,targetAx.Position,dcm_obj,targetCoords,X,labs));
else
    set(fig,'units','normalized','WindowButtonDownFcn',{@circIfIn,targetAx.Position}, ...
        'WindowButtonUpFcn',@circOff)
    lh = addlistener(fig,'CurrentPoint','PostSet', ...
        @(src,evnt)updateTextCnt(src,evnt,dataPlot,targetAx,targetCoords,X,labs));
end


%% callbacks
% for the cursor case
    function updateTextCsr(~,~,varargin) 
        datPlt = varargin{1};
        axLims = varargin{2};
        dcm = varargin{3};
        targDat = varargin{4};
        givenDat = varargin{5};
        labels = varargin{6};
        dcm_inf = getCursorInfo(dcm);
        if isWithin(fig.CurrentPoint,axLims) && ~isempty(dcm_inf)
            cp = dcm_inf.Position;
            
            dd = [targDat - cp(1,1:size(targDat,2))]; % assume we want the point closest to the screen
            dist2pnt = sqrt(diag(dd*dd'));
            [~,ind] = min(dist2pnt); % will this always be unique???
            if isa(datPlt,'matlab.graphics.primitive.Image')
                datPlt.CData = givenDat{ind};
            else
                datPlt.YData = givenDat(:,ind)';
            end
            title(get(datPlt,'Parent'),labels{ind})
        end
    end

% for the continuous case
    function updateTextCnt(~,~,varargin)
        datPlt = varargin{1};
        targAx = varargin{2};
        targDat = varargin{3};
        givenDat = varargin{4};
        labels = varargin{5};
        if isWithin(fig.CurrentPoint,targAx.Position)
            cp = targAx.CurrentPoint;
            
            dd = [targDat - cp(1,1:size(targDat,2))]; % assume we want the point closest to the screen
            dist2pnt = sqrt(diag(dd*dd'));
            [~,ind] = min(dist2pnt); % will this always be unique???
            if isa(datPlt,'matlab.graphics.primitive.Image')
                datPlt.CData = givenDat{ind};
            else
                datPlt.YData = givenDat(:,ind)';
            end
            title(get(datPlt,'Parent'),labels{ind})
        end
    end

    function circIfIn(src,~,targLims) % button down function
        if isWithin(src.CurrentPoint,targLims)
            src.Pointer = 'circle';
            src.WindowButtonMotionFcn = {@takeBreath,targLims};
        end
    end

    function circOff(src,~) % button up function
        src.Pointer = 'arrow';
        src.WindowButtonMotionFcn = '';
    end

% general
    function takeBreath(src,~,targLims) % motion function
        if isWithin(src.CurrentPoint,targLims)
            pause(0.01); % so that we don't try to update the plots too fast
        end
    end

end

%% helpers
function isIn = isWithin(pnt,bnds)
% decide whether point in figure coordinates is within given axis position.
% where axis position is [left bottom width height] in figure units

isIn = all(pnt > bnds([1 2])) && all(pnt < (bnds(3:4)+bnds(1:2)));

end

function type = whichClass(thisObj)
% strange workaround to determine class of a gaphics object -- is there a
% native function for this? I don't think so, and this is not super robust
% at all so do replace it if you find one later.

% use isa('matlab.graphics.chart.primitive.Scatter')

type = cell(1,length(thisObj));
for ii = 1:length(thisObj)
    try
        get(thisObj(ii),'fakeproperty');
    catch err
        [~,tmp1] = regexp(err.message,'property on the ');
        tmp2 = regexp(err.message,'class');
        type{ii} = err.message(tmp1+1:tmp2-2);
    end
end

if length(type) == 1, type = type{:}; end

end

function tabulaRasa(f)

try
    predecessors = f.AutoListeners__;
    for ii = 1:length(predecessors)
        delete(predecessors{ii});
    end
catch
end

end