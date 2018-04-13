function plot2DTuning(ba,whichVar,bins,pp)
% plot2DTuning(BA,whichVar,bins[,params])
%
% Takes a binned array of spike counts and makes a nice representation of
% the cell's tuning based on those counts. See 'psthAndBa' function, which
% can generate the first input.
%
% Inputs
%   - BA [nTrials,nBins]: Binned spike counts for each trial
%   - whichVar [nTrials,2]: Event labels for each trial, where the first
%   column is for one variable and the second for another (e.g., first is
%   the strength of one stimulus, second is strength of another)
%   - bins [nBin,1]: time label for each bin
%   - params (struct): Options with fields:
%       > majorVar (1 or 2): which variable should go on the y-axis, and be
%       the outer variable on the PSTH plots. Number is column of whichVar.
%           Default: 1
%       > varlabs {var1,var2}: Labels for each column in whichVar
%           Default: {'Var1,'Var2'}
%       > respWin [nBin,1] logical: Which bins to include in 'response'
%       calculation (e.g. (bins <= 0) for a causal index).
%           Default: (bins >=0)
%       > useDiff (logical): Whether to subtract baseline in calculation of
%       'response'. Baseline here means the counts that aren't in the
%       response window.
%           Default: true
%       > whichAx (5x1 axes): Supply the axis handles for each of the
%       output plots. The order is:
%           [tuningcurve,topMarginal,rightMarginal,meanPSTH,allPSTH]
%           By default it will make these axes on the current figure.
%
% M.A., UCL 2018

% defaults
whichAx = [];
majorVar = 1;
varlabs = {'Var 1' 'Var 2'};
respWin = (bins>=0); % when to compute 'responsiveness'
useDiff = true;

if exist('pp','var') % replace by input parameters
    if isfield(pp,'whichAx'), whichAx = pp.whichAx; end
    if isfield(pp,'majorVar'), majorVar = pp.majorVar; end
    if isfield(pp,'varlabs'), varlabs = pp.varlabs; end
    if isfield(pp,'respWin'), respWin = pp.respWin; end
    if isfield(pp,'useDiff'), useDiff = pp.useDiff; end
end

minorVar = setdiff([1 2],majorVar);
win = [min(bins) max(bins)];
binsize = mean(diff(bins));
dt_resp = range(bins(respWin));

if size(whichVar,1) ~= size(ba,1)
    whichVar = whichVar';
end
var1 = whichVar(:,majorVar); % the variable on the y axis
var1Lab = varlabs{majorVar};
var2 = whichVar(:,minorVar); 
var2Lab = varlabs{minorVar};

% compute responses
v1 = unique(var1); 
v2 = unique(var2);
nVar1 = length(v1);
nVar2 = length(v2);
resp = zeros(nVar1,nVar2);
errs = zeros(nVar1,nVar2);
psth = zeros(nVar2*nVar1,length(bins));
nt = zeros(nVar1,nVar2);
dingle = zeros(nVar2*nVar1,2); % for debugging
dongle = zeros(nVar1,nVar2,2);
for ii = 1:nVar1
    for jj = 1:nVar2
        theseEvents = (var1 == v1(ii)) & (var2 == v2(jj));
        p = mean(ba(theseEvents,:)./binsize,1);
        
        if useDiff
            r0 = mean(sum(ba(theseEvents,~respWin),2)/dt_resp);
        else
            r0 = 0;
        end
        spkCounts = sum(ba(theseEvents,respWin),2)/dt_resp;
        
        nt(ii,jj) = sum(theseEvents);
        resp(ii,jj) = mean(spkCounts - r0);
        errs(ii,jj) = std(spkCounts - r0)./sqrt(sum(theseEvents));
        psth(jj + nVar2*(ii-1),:) = p;
        dingle(jj + nVar2*(ii-1),:) = [v1(ii) v2(jj)];
        dongle(ii,jj,1) = v1(ii);
        dongle(ii,jj,2) = v2(jj);
    end
end

%% plot
gcf;
if isempty(whichAx)
    surfAx = subtightplot(3,6,[7 8 13 14],[0.05 0.01],[0.1 0.05],[0.07 0.1]);
    marg2 = subtightplot(3,6,[1 2],[0.01 0.01],[0.05 0.05],[0.07 0.1]);
    marg1 = subtightplot(3,6,[9 15],[0.05 0.0005],[0.1 0.05],[0.05 0.05]);
    
    tcAx = subtightplot(3,6,[4 5 6],[0.01 0.1]);
    psthAx = subtightplot(3,6,[10 11 12  16 17 18],[0.05 0.1],[0.1 0.05],[0.05 0.05]);
else
    surfAx = whichAx(1);
    marg2 = whichAx(2);
    marg1 = whichAx(3);
    
    tcAx = whichAx(4);
    psthAx = whichAx(5);
end


% tuning surface with marginals
[~,p1] = sort(v1(:),1,'descend'); % plotting order
[~,p2] = sort(v2(:));
    % full tuning surface

imagesc_but_good(surfAx,v2,v1(p1),resp(p1,p2))
hold(surfAx,'on');
tx = strcat('n= ',(split(num2str(reshape(nt(p1,p2),nVar1*nVar2,1)'))));
text(surfAx,repelem(1:nVar2,1,nVar1),repmat(1:nVar1,1,nVar2),tx,'horizontalalignment','center')
xlabel(surfAx,var2Lab); ylabel(surfAx,var1Lab)
hold(surfAx,'off');
    % upper marginal
errorbar(marg2,nanmean(resp(p1,p2),1),nanmean(errs(p1,p2),1),'-o')
xticks(marg2,sort(p2));
xticklabels(marg2,[])
xlim(marg2,[0.5,nVar2+0.5])
ylabel(marg2,'Mean FR (Hz)')
    % right marginal
errorbar(marg1,nanmean(resp(p1,p2),2),(p1),[],[],...
    nanmean(errs(p1,p2),2),nanmean(errs(p1,p2),2),'-o')
yticks(marg1,sort(p1))
yticklabels(marg1,[])
ylim(marg1,[0.5,nVar1+0.5])
xlabel(marg1,'Mean FR (Hz)')
    % make axes equal
marglims = [min([ylim(marg2) xlim(marg1)]),max([ylim(marg2) xlim(marg1)])];
marg2.YLim = marglims;
marg1.XLim = marglims;

% response timecourse
B = binmat(size(psth,2),10)/10; % bin the psth
[~,p2] = sort(v2(:),1,'descend');
psthOrder = repmat(p2,nVar1,1) + nVar2*(repelem(p1,nVar2)-1);
    % mean psth
plot(tcAx,bins*B',nanmean(psth(:,:),1)*B')
hold(tcAx,'on')
plot(tcAx,[0 0],tcAx.YLim,'--','color','r')
plot(tcAx,bins(respWin),max(ylim(tcAx))*ones(sum(respWin),1),'r','linewidth',4)
text(tcAx,mean(bins(respWin)),max(ylim(tcAx)),'''Response''','fontsize',12,...
    'verticalalignment','top','horizontalalignment','center')
xticks(tcAx,[])
ylabel(tcAx,'Firing Rate (Hz)')
hold(tcAx,'off')
    % all psth
imagesc(psthAx,bins,[],psth(psthOrder,:)*B')
psthAx.YTick = [2:nVar1:nVar1*nVar2] + 0.5;
hold(psthAx,'on');
line_y = repmat([nVar2 nVar2 NaN],1,nVar1) + nVar2*repelem(0:(nVar1-1),3) + 0.5;
plot(psthAx,repmat([xlim NaN],nVar1),line_y,'--','color','w','linewidth',2)
plot(psthAx,[0 0],psthAx.YLim,'--','color','r','linewidth',2)
yticklabels(psthAx,v1(p1))
yy1 = ylabel(psthAx,var1Lab);
set(yy1,'Units','Normalized','Position',[-0.1, 0.5, 0])
xlabel(psthAx,'time')
hold(psthAx,'off');
box(psthAx,'off');
foo = axes('Position',get(psthAx,'Position'),'color','none');
set(foo,'yaxislocation','right')
foo.YLim = [0 nVar1*nVar2];
foo.YTick = repmat([0 (nVar2-1)],1,nVar1) + nVar2*repelem(0:(nVar1-1),2) + 0.5;
foo.YTickLabels = repmat([min(v2) max(v2)],1,nVar1);
foo.XTick = [];
yy = ylabel(foo,var2Lab,'rotation',270);
set(yy,'Units','Normalized','Position',[1.1, 0.5, 0]);

linkaxes([psthAx,tcAx],'x')


