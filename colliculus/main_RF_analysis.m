muafeats = loadVar('C:\DATA\Spikes\SC_popstructs\MUA_ephys_features.mat','feats');
allfeats = loadVar('C:\DATA\Spikes\SC_popstructs\all_ephys_features.mat','feats');
clu_fits = loadVar('C:\DATA\Spikes\SC_popstructs\all_RF_fits.mat','clu_fits');
% q_fits = loadVar('C:\DATA\Spikes\SC_popstructs\RF_fits_Q.mat','clu_fits');

b = length(muafeats.datasets) + 1;
wc2id = @(WC) WC*[b^2 b^1 b^0]'; % from [CID,tag,dset] -> unique CID
id2wc = @(id) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)]; % inverse

% ----------------------
% fix situation where the two dataset structures are dissimilar
[~,g1,g2] = dbComp(muafeats.datasets,clu_fits.datasets);
tmp_cfwc = clu_fits.whichCell;
tmp_muawc = muafeats.whichCell;
[foo1,locb] = ismember(tmp_muawc(:,3),g1);
tmp_muawc(:,3) = locb;
% ----------------------

cfid = wc2id(tmp_cfwc);
muaid = wc2id(tmp_muawc);

%% parameters
nongaus_thresh = 0.005;
bilat_thresh = 0.013;
goodFit_thresh = 0.03;
sigma_thresh = 7; % minimum z-score of RF peak

dontUse = {'fw3m','capacitative','acg half-width','acg 2/3-width','mean rate','spontaneous rate','hds isi'};
goodfeats = ~contains(lower(muafeats.labels),dontUse); % some features are fishy/not useful
goodclu = ~any(isnan(muafeats.features(:,goodfeats)),2); % some cells have NaN in one or more features, ignore them
%% make cluster groups
R2 = mean(clu_fits.score,3);
isModulated = [clu_fits.rfStats(:).z_peak] > sigma_thresh;
isFitWell = clu_fits.empiricalScore > 0;
chosenCells = isModulated(:) & isFitWell(:);

inBoth = ismember(muaid,cfid(chosenCells));
inSCs = muaid(inBoth & (muafeats.depth(:,3) >= 0));
inSCi = muaid(inBoth & (muafeats.depth(:,3) < 0) & (muafeats.depth(:,3) >= -300));
inSCd = muaid(inBoth & (muafeats.depth(:,3) < -300));

sens = ismember(cfid,inSCs);
mot = ismember(cfid,inSCi);
othr = ismember(cfid,inSCd);

bilatBest = (R2(:,4) - max(R2(:,1:3),[],2) > bilat_thresh)  & (chosenCells );
% bilatBest = (R2(:,4) - max(R2(:,1:3),[],2) > r2_thresh) & (isModulated & isFitWell);
unilatBest = (~bilatBest) & (chosenCells );
nongausBest = (max(R2(:,2:3),[],2) - R2(:,1) > nongaus_thresh) & (chosenCells );
gausBest = (~nongausBest) & (chosenCells );

nClu = sum(chosenCells );
%% get some RF parameters
std_main = nan(length(cfid),2);
mainCenter = nan(length(cfid),2);
cntrs = nan(length(cfid),4);
for iClu = cfid(: )'
    std_main(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{1}([3 5]);
    mainCenter(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{1}([2 4]);
    cntrs(cfid == iClu,1:2) = clu_fits.bestParams{cfid == iClu}{1}([2 4]);
end

std_anti = nan(length(cfid),2);
antiCenter = nan(length(cfid),2);
antiAmps = nan(length(cfid),1);
amps = nan(length(cfid),2);
cntrs_distance = nan(length(cfid),1);
for iClu = cfid(:)'
    if bilatBest(cfid == iClu)
        std_anti(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{2}{2}([3 5]);
        antiCenter(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{2}{2}([2 4]);
        antiAmps(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{2}{2}(1);
    end
    if nongausBest(cfid == iClu)
        amps(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{1}([1 6]);
        if length(clu_fits.bestParams{cfid == iClu}{1}) > 9
            cntrs(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{1}([2 4 9 10]);
            dd = clu_fits.bestParams{cfid == iClu}{1}([9 10])' - mainCenter(cfid == iClu,:);
            cntrs_distance(cfid == iClu) = norm(dd);
        else
            cntrs(cfid == iClu,:) = clu_fits.bestParams{cfid == iClu}{1}([2 4 2 4]);
            cntrs_distance(cfid == iClu) = 0;
        end
    end
end

nongausBest = nongausBest & (cntrs_distance(:) <= mean(std_main,2));
gausBest = (~nongausBest) & (chosenCells );

%% distribution of various things, separated by groups
% grplabs = {'sSC','iSC','dSC and below'};
grplabs = {'Single gaussian','Center-suppressive'};
% grplabs = {'Unilateral','Bilateral'};
% grps = [chosenCells ~chosenCells];
% grps = [sens mot othr];
% grps = [gausBest nongausBest];
grps = [unilatBest bilatBest] ;
grps = grps(ismember(cfid,muaid),:);
% quantity = (R2(:,4)- max(R2(:,1:3),[],2));
% quantity = (max(R2(:,2:3),[],2) - R2(:,1));
% quantity = max(R2(:,:),[],2) - clu_fits.rfR2;
% quantity = nanmean(std_main,2);
% quantity = clu_fits.clusterInfo(:,3);
quantity = muafeats.depth(ismember(muaid,cfid),3);
foo = logical(sum(grps,2));
e = linspace(quantile(quantity(foo),0.0001),quantile(quantity(foo),0.99),34);

subplot(2,1,2); hold on
alln = zeros(size(grps,2),length(e)-1);
kyu = {};
for ii = 1:size(grps,2)
    dis = grps(:,ii);
    kyu{ii} =quantity(dis);
     histogram(quantity(dis),e,'normalization','probability','edgecolor','none'); % ,'DisplayStyle','stairs','linewidth',2
%     n = histcounts(quantity(dis),e);
%     alln(ii,:) = n/sum(dis);
%     [n,e] = ecdf(quantity(dis));
%     plot(e,n,'linewidth',2)
end

% figure;
% barh(e(1:end-1),alln(1,:),1,'edgecolor','none','facealpha',0.5)
% hold on
% barh(e(1:end-1),alln(2,:),1,'edgecolor','none','facealpha',0.5)
% bar(repmat(e(1:end-1),size(grps,2),1)',alln',1,'stacked','edgecolor','none')
% colormap flag
% xlabel('Depth from stratum opticum')
% ylabel('Fraction of visual cells')

legend(grplabs)

% bar(e(1:end-1),alln(1,:)./nCluInBin,1,'edgecolor','none','facecolor','b','facealpha',0.7)
% hold on
% bar(e(1:end-1),alln(2,:)./nCluInBin,1,'edgecolor','none','facecolor','r','facealpha',0.7)
% bar(e(1:end-1),alln(3,:),1,'edgecolor','none','facecolor','y','facealpha',0.7)

%% PCA (and plots)
trueMUA = ~ismember(muafeats.whichCell,allfeats.whichCell,'rows');

goodclu_all = ~any(isnan(allfeats.features(:,goodfeats)),2);
X = allfeats.features(goodclu_all,goodfeats);
zX = (X-nanmean(X,1));
zX = bsxfun(@rdivide,zX,nanstd(X,1));

[W_all, P, V] = pca(zX);
V = V./sum(V);

X_mua = muafeats.features(trueMUA & goodclu,goodfeats);
zX_mua = (X_mua-nanmean(X_mua,1));
zX_mua = bsxfun(@rdivide,zX_mua,nanstd(X_mua,1));

P_mua = zX_mua*W_all;
V_mua = diag(cov(P_mua))./sum(diag(cov(zX_mua)));

figure;plot(V,'k','linewidth',2)
hold on
plot(V_mua,'--','color','k','linewidth',2)
legend({'Single units','MUA'})
xlabel('PC','fontsize',13)
ylabel('% var','fontsize',14)
xlim([1 length(V)])
box off

figure;
imagesc(W_all(:,1:3))
colormap(redblue)
yticks(1:sum(goodfeats))
yticklabels(muafeats.labels(goodfeats))
c = colorbar;
c.Label.String = 'weight';
c.Label.FontSize = 15;
c.FontSize = 15;
caxis([-0.5 0.5]);
xlabel('PC')

figure;
scatter3(P_mua(:,1),P_mua(:,2),P_mua(:,3),36,muafeats.depth(goodclu &trueMUA,3),'.');
xlabel(['PC1 (' num2str(100*V(1),2) '% var)'],'fontsize',14)
ylabel(['PC2 (' num2str(100*V(2),2) '% var)'],'fontsize',14)
zlabel(['PC3 (' num2str(100*V(3),2) '% var)'],'fontsize',14)
xlabel(['PC1'],'fontsize',14)
ylabel(['PC2'],'fontsize',14)
zlabel(['PC3'],'fontsize',14)
c = colorbar;
c.Label.String = 'Depth from S.O.';
c.FontSize = 12;
caxis([-2000 900]) % make SC bottom around the yellow/blue border
view([-125.9 54])
grid off

% RF groups in this space
grp1 = ismember(muaid(trueMUA & goodclu),cfid(chosenCells));
% grp2 = ismember(muaid(trueMUA & goodclu),cfid(~chosenCells));
grp2 = false(size(grp1));

figure;
scatter3(P_mua(grp1,1),P_mua(grp1,2),P_mua(grp1,3),14,'k','filled')
hold on
% scatter3(P_mua(grp2,1),P_mua(grp2,2),P_mua(grp2,3),14,'b','filled')
scatter3(P_mua(~(grp1 | grp2),1),P_mua(~(grp1 | grp2),2),P_mua(~(grp1 | grp2),3),24, ...
    [0.5 0.5 0.5],'.');

xlabel(['PC1 (' num2str(100*V(1),2) '% var)'])
ylabel(['PC2 (' num2str(100*V(2),2) '% var)'])
zlabel(['PC3 (' num2str(100*V(3),2) '% var)'])
legend({'Nonvisual','Visual'})
grid off

%% LDA (and plots)
useThese = ismember(muaid,cfid) & goodclu & inBoth;
grplabs = {'Single gaussian','Center-suppressive'};
% grplabs = {'Unilateral','Bilateral'};
% grp1 = ismember(muaid,cfid(gausBest));
% grp2 = ismember(muaid,cfid(nongausBest));
grp1 = ismember(muaid,cfid(unilatBest));
grp2 = ismember(muaid,cfid(bilatBest));
% grp1 = ismember(muaid,cfid(chosenCells));
% grp2 = ismember(muaid,cfid(~chosenCells));
grp1 = grp1(useThese);
grp2 = grp2(useThese);

X = muafeats.features(useThese,goodfeats);
zX = (X-nanmean(X,1));
zX = bsxfun(@rdivide,zX,nanstd(X,1)); % standardize for the population

X_vis = zX(:,:);
C = nan(size(X_vis,1),1); % what to discriminate by
C(grp1) = 1;
C(grp2) = 2;

W = goodLDA(X_vis,C);
L = X_vis*W';
%%
bns = linspace(quantile(L,0.001),quantile(L,0.98),11);
figure;histogram(L(grp1),bns,'normalization','probability','edgecolor','none')
hold on
histogram(L(grp2),bns,'normalization','probability','edgecolor','none')
xlabel('Discriminant value')
ylabel('Fraction of RF type')
legend(grplabs)

figure;barh(W)
yticklabels(muafeats.labels(goodfeats))
set(gca,'ydir','reverse')
set(gca,'YAxisLocation','right')
box off
xlabel('Weight')
xlim([-1 1])

% maybe see it in 2d?
otherdims = null(W);
[~,pea] = pca(otherdims); % most variable directions ortho to LD1
L_2d = X_vis*pea(:,1);

%% look at features interactively (got to mess with zoom first)
grp1 = ismember(muaid,cfid(nongausBest)) & goodclu;% & ismember(muaid,wc2id(allfeats.whichCell));
grp2 = ismember(muaid,cfid(gausBest))& goodclu;% & ismember(muaid,wc2id(allfeats.whichCell));

f2p = [9,14,1]; % features to plot

% theseACG = [flipud(muafeats.ACG(grp1,:)'); muafeats.ACG(grp1,2:end)';];
theseACG = [flipud(muafeats.ACG(grp2,:)'); muafeats.ACG(grp2,2:end)';];
 
theseBins = [-fliplr(muafeats.ACG_bins) muafeats.ACG_bins(2:end)];

figure('units','normalized','position',[0.1182 0.4361 0.6505 0.4259]);
ax1 = subplot(1,2,1);
scatter3(muafeats.features(grp1,f2p(1)),muafeats.features(grp1,f2p(2)),muafeats.features(grp1,f2p(3)),'filled')
hold on
scatter3(muafeats.features(grp2,f2p(1)),muafeats.features(grp2,f2p(2)),muafeats.features(grp2,f2p(3)),14,[0.5 0.5 0.5],'filled')
xlabel([muafeats.labels(f2p(1)) ' (sec)'])
ylabel([muafeats.labels(f2p(2)) ' (sec)'])
zlabel([muafeats.labels(f2p(3)) ' (sec)'])
view([60,13.2])
legend({'Bullseye','Simple'},'location','best')

subplot(1,2,2)
% badplot(ax1,theseACG,theseBins,muaid(grp1),true)
badplot(ax1,theseACG,theseBins,muaid(grp2),true)
xlabel('Lag (sec)')
ylabel('Conditional intensity (sp/sec)')
xlim([-0.12 0.12])

disp('Use datacursor to select cells')

%% look at feature for single cell
thisclu = 279608;
thislag = muafeats.features(muaid == thisclu,9);

figure
ax = subplot(2,1,1);
bar(theseBins',theseACG(:,muaid(grp1) == thisclu),1,'facecolor',[0.3 0.3 0.3])
xlim([-0.12 0.12])
xlabel('Lag (sec)','fontsize',13)
ylabel('Rate (sp/sec)','fontsize',13)
box off

mlind = find(theseBins >= thislag,1,'first');
ymin = theseACG(mlind,muaid(grp1) == thisclu);
ymax = ymin + 5;

% annotations are annoyingly in 'normalized' coordinates, so I need to
% convert from actual coordinates to values normalized for the figure, i.e.
% (0,0) is lower left, (1,1) upper right. 'pos' tells you the position of
% the axes in these 'normalized' coordinates.
pos = ax.Position;
x2norm = @(x)(x-ax.XLim(1))/diff(ax.XLim)*pos(3) + pos(1);
y2norm = @(y)(y-ax.YLim(1))/diff(ax.YLim)*pos(4) + pos(2);
annotation('textarrow',x2norm([thislag thislag]),y2norm([ymax ymin]),'String',[num2str(100*thislag,3) 'ms'])

% spike shape
foo = id2wc(thisclu);
CID = foo(1);
k = foo(3); t = foo(2);
if (k~= prevk) || (t~=prevt)
    [spks,dataDir,ksDir] =  quickLoadData(muafeats.datasets,foo,'spikes');
    [~,max_site_temps] = max(max(abs(spks.tempsUnW),[],2),[],3); % not zero-indexed
    max_site = max_site_temps(spks.cluTemps(spks.cids+1)+1); % need to add 1 for indexing, both times
    prevk = k; prevt = t;
end
gwfparams.dataDir = dataDir; % for getting mean waveforms
gwfparams.ksDir   = ksDir;     % takes a loooong time
gwfparams.fileName = spks.dat_path;
gwfparams.dataType = spks.dtype;
gwfparams.nCh = spks.n_channels_dat;
gwfparams.wfWin = [-31 35];
gwfparams.nWf = 300;
gwfparams.spikeTimes = round(spks.st(spks.clu == CID)*spks.sample_rate);
gwfparams.spikeClusters = spks.clu(spks.clu == CID);
gwfparams.nMaxCh = [1 1];
gwfparams.maxSite = max_site(spks.cids == CID);
wfMainChan = getWaveFormsMA(gwfparams);
wf = spks.gain*permute(wfMainChan.waveForms,[4 2 3 1]);

subplot(2,1,2)
plot(wf,'color',[0.5 0.5 0.5])
axis off

%% ???
% cfidepths = zeros(sum(deez),1);
% j = 1;
localRoot = 'C:\DATA\Spikes\';
prevk = NaN; prevt = NaN;
fooRF = cell(sum(chosenCells),1);
fooSts = struct('timeBins',{},'peakZscore',{},'scalar',{},'bsl',{},'p_zscore',{},'metaZ',{});
theseid = cfid(chosenCells);
for iClu = 1:sum(chosenCells)
    foo = id2wc(theseid(iClu));
    k = foo(3); t = foo(2);
    CID = foo(1);
    
    if (k~= prevk) || (t~=prevt)
        thisTag = clu_fits.datasets(k).tags{t};
        
        thisExp =  [db(k).mouse_name '_' db(k).date '_' thisTag];
        %         if verbose, dispstat(['Starting ' thisExp],'keepthis','timestamp'); end
        
        % load everything for that probe
        [dsetFolders,dataDir] = expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [dataDir '\sorting\'];
        
        spks = loadKSdir(ksDir);
        
        snname = [localRoot dsetFolders '\sparse_noise_RFs.mat'];
        snrf = loadVar(snname,'snrf');
        noiseStart = min(snrf.stimTimes_local);
        noiseEnd = max(snrf.stimTimes_local);
        
        prevk = k; prevt = t;
    end
    
    deezspks = (spks.clu == CID) & (spks.st>noiseStart & spks.st<noiseEnd);
    p = struct('useSVD',true,'makePlots',false,'fit2Dgauss',false,'nShuffle',[]);
    p.binSize = 0.01;
    [RF,sts] = sparseNoiseRF_MA(spks.st(deezspks), snrf.stimTimes_local,snrf.stimPosition,p);
    fooRF{iClu} = RF;
    fooSts(iClu) = sts;
end

%%
thisclu = 16932;
foo = id2wc(thisclu);
k = foo(3); t = foo(2);
CID = foo(1);

if (k~= prevk) || (t~=prevt)
    thisTag = clu_fits.datasets(k).tags{t};
    
    thisExp =  [db(k).mouse_name '_' db(k).date '_' thisTag];
    %         if verbose, dispstat(['Starting ' thisExp],'keepthis','timestamp'); end
    
    % load everything for that probe
    [dsetFolders,dataDir] = expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
    ksDir = [dataDir '\sorting\'];
    
    spks = loadKSdir(ksDir);
    
    snname = [localRoot dsetFolders '\sparse_noise_RFs.mat'];
    snrf = loadVar(snname,'snrf');
    noiseStart = min(snrf.stimTimes_local);
    noiseEnd = max(snrf.stimTimes_local);
    
    prevk = k; prevt = t;
end
deezspks = (spks.clu == CID) & (spks.st>noiseStart & spks.st<noiseEnd);
p = struct('useSVD',true,'makePlots',false,'fit2Dgauss',false,'nShuffle',[]);
p.binSize = 0.01;
[RF,sts,PSTH] = sparseNoiseRF_MA(spks.st(deezspks), snrf.stimTimes_local,snrf.stimPosition,p);
[allpos,~,pos_grps] = unique(snrf.stimPosition,'rows');

figure
imagesc(snrf.XPos,-snrf.YPos,RF)
q = abs(quantile(abs(RF(:)),0.98));
caxis([-q q])

%% x-corrs
bs = 0.01; % bin size for psth
vis_win = [-0.05 0.2]; % window for each kind of event

ccg_bs = 0.01;
ccg_win = 0.1; % +/- from 0, in sec
min_spont = 60; % minimum interval to be considered 'spontaneous', in sec



%%
% [X,Y] = meshgrid(snrf.XPos,snrf.YPos);
dis = clu_fits.bestParams{cfid == thisclu}{1};
mu = cntrs(cfid == thisclu,:);
M = [cos(dis(end)) -sin(dis(end)); sin(dis(end)) cos(dis(end))];

foo_main = M*[real(1.5*dis(3)*exp(1i*[-1:0.1:2*pi])); ...
                imag(1.5*dis(5)*exp(1i*[-1:0.1:2*pi]))] + [mu(1);-mu(2)];
cir = foo_main;
if nongausBest(cfid == thisclu)
    foo_cntr =  M*[real(1.5*dis(7)*exp(1i*[-1:0.1:2*pi])); ...
        imag(1.5*dis(8)*exp(1i*[-1:0.1:2*pi]))] + [mu(3);-mu(4)];
    cir = [cir nan(2,1) foo_cntr];
end
           

if bilatBest(cfid == thisclu)
    dis = clu_fits.bestParams{cfid == thisclu}{2}{2};
    foo_anti = M*[real(1.5*dis(3)*exp(1i*[-1:0.1:2*pi])); ...
                imag(1.5*dis(5)*exp(1i*[-1:0.1:2*pi]))] + [dis(2);-dis(4)];
    cir = [cir nan(2,1) foo_anti];
end

hold on
plot(cir(1,:),cir(2,:),'--k','linewidth',2)
%%
figure('Position',[27 520 524 265])
% subplot(2,2,1)
% imagesc(snrf.XPos,-snrf.YPos,RF)
% pause
loc1 = [-84, -25];
% loc2 = [114 -11];
% loc3 = [-114 4];
% title(['Clu: ' num2str(foo)])

subtightplot(4,1,2:4)
dis = find(round(allpos(:,2)) == loc1(1) & round(allpos(:,1)) == loc1(2));
[psth,bins,rasx,rasy] = psthRasterAndCounts(spks.st(deezspks),snrf.stimTimes_local(pos_grps == dis),[-0.05 0.2],0.01);
plot(rasx,rasy,'k','linewidth',2)
xlim([-0.05 0.2])
% title(['Location: ' num2str(loc1)])
hold on
plot([0 0],ylim,'-r')
plot([0.0156 0.0156],ylim,'-b')
axis off

subtightplot(4,1,1,[0.001, 0.01])
plot(bins,psth,'k','linewidth',2)
hold on
plot([0 0],ylim,'-r')
plot([0.0156 0.0156],ylim,'-b')
midpnt = mean(ylim);
plot([min(xlim) min(xlim)],[midpnt-5 midpnt+5],'k','linewidth',2)
xlim([-0.05 0.2])
axis off

% figure;
% plot(sts.timeBins,sts.timeCourse,'k','linewidth',2)

% subplot(2,2,3)
% dis = find(round(allpos(:,2)) == loc2(1) & round(allpos(:,1)) == loc2(2));
% [~,~,rasx,rasy] = psthRasterAndCounts(spks.st(deezspks),snrf.stimTimes_local(pos_grps == dis),[-0.05 0.2],0.01);
% plot(rasx,rasy,'k')
% title(['Location: ' num2str(loc2)])
% hold on
% plot([0 0],ylim,'-r')
% plot([0.0156 0.0156],ylim,'-b')
% 
% subplot(2,2,4)
% dis = find(round(allpos(:,2)) == loc3(1) & round(allpos(:,1)) == loc3(2));
% [~,~,rasx,rasy] = psthRasterAndCounts(spks.st(deezspks),snrf.stimTimes_local(pos_grps == dis),[-0.05 0.2],0.01);
% plot(rasx,rasy,'k')
% title(['Location: ' num2str(loc3)])
% hold on
% plot([0 0],ylim,'-r')
% plot([0.0156 0.0156],ylim,'-b')


%%
% [~,inds] = sort(R2(bilatBest,4) - max(R2(bilatBest,1:3),[],2),'descend');
[~,inds] = sort(clu_fits.clusterInfo(bilatBest,3),'descend');
bills = cfid(bilatBest);
bills = bills(inds(1:30));
bills = reshape(reshape(bills,10,3)',30,1);
figure
for ii = 1:30
    subtightplot(10,3,ii,[0 0],[0 0],[0 0]);
    imagesc(clu_fits.RF{cfid == bills(ii)});
    q = abs(quantile(abs(clu_fits.RF{cfid == bills(ii)}(:)),0.97));
    caxis([-q q])
    
    hold on
    rfx = size(clu_fits.RF{cfid == bills(ii)},2) + 0.5;
    rfy = size(clu_fits.RF{cfid == bills(ii)},1) + 0.5;
    xmid = size(clu_fits.RF{cfid == bills(ii)},2)/2 + 0.5;
    if clu_fits.clusterInfo(cfid == bills(ii),3) > 0
        plot([0.5 rfx rfx 0.5 0.5],[repelem([0.5 rfy],2) 0.5],'color','b')
    else
        plot([0.5 rfx rfx 0.5 0.5],[repelem([0.5 rfy],2) 0.5],'color','r')
    end
    plot([xmid xmid],[0.5 rfy],'-k')
    
%     dis = clu_fits.bestParams{cfid ==  bills(ii)}{1};
%     mu = cntrs(cfid ==  bills(ii),:);
%     M = [cos(dis(end)) -sin(dis(end)); sin(dis(end)) cos(dis(end))];
%     
%     foo_main = M*[real(1.5*dis(3)*exp(1i*[-1:0.1:2*pi])); ...
%         imag(1.5*dis(5)*exp(1i*[-1:0.1:2*pi]))] + [mu(1);-mu(2)];
%     cir = foo_main;
%     if nongausBest(cfid ==  bills(ii))
%         foo_cntr =  M*[real(1.5*dis(7)*exp(1i*[-1:0.1:2*pi])); ...
%             imag(1.5*dis(8)*exp(1i*[-1:0.1:2*pi]))] + [mu(3);-mu(4)];
%         cir = [cir nan(2,1) foo_cntr];
%     end
%     
%     if bilatBest(cfid ==  bills(ii))
%         dis2 = clu_fits.bestParams{cfid ==  bills(ii)}{2}{2};
%         foo_anti = M*[real(1.5*dis2(3)*exp(1i*[-1:0.1:2*pi])); ...
%             imag(1.5*dis2(5)*exp(1i*[-1:0.1:2*pi]))] + [dis2(2);-dis2(4)];
%         cir = [cir nan(2,1) foo_anti];
%     end
% 
%     hold on
%     plot(cir(1,:),cir(2,:),'--k','linewidth',2)
    axis equal;axis off;
end

%%
% [~,inds] = sort(max(R2(nongausBest,2:3),[],2) - R2(nongausBest,1),'descend');
[~,inds] = sort(clu_fits.clusterInfo(nongausBest,3),'descend');
dougs = cfid(nongausBest);
dougs = dougs(inds(1:50));
dougs = reshape(reshape(dougs,10,5)',50,1);
figure
for ii = 1:50
    subtightplot(10,5,ii,[0 0],[0 0],[0 0]);
    imagesc(clu_fits.RF{cfid == dougs(ii)});
    q = abs(quantile(abs(clu_fits.RF{cfid == dougs(ii)}(:)),0.97));
    caxis([-q q])
    axis equal;axis off;
    hold on
    rfx = size(clu_fits.RF{cfid == dougs(ii)},2) + 0.5;
    rfy = size(clu_fits.RF{cfid == dougs(ii)},1) + 0.5;
    xmid = size(clu_fits.RF{cfid == dougs(ii)},2)/2 + 0.5;
    if clu_fits.clusterInfo(cfid == dougs(ii),3) > 0
        plot([0.5 rfx rfx 0.5 0.5],[repelem([0.5 rfy],2) 0.5],'color','b')
    else
        plot([0.5 rfx rfx 0.5 0.5],[repelem([0.5 rfy],2) 0.5],'color','r')
    end
    plot([0.5 rfx rfx 0.5 0.5],[repelem([0.5 rfy],2) 0.5],'-k')
    plot([rfx/2 rfx/2],[0.5 rfy],'-k')
end

%% location of training stimuli
db = clu_fits.datasets;
numdates = 7;

bills = clu_fits.whichCell(bilatBest,3);
[thesemice,inds,kays] = unique({db(bills).mouse_name});
datalocs = {db(bills).dataServer};
datalocs = datalocs(inds);
cwexp = {db(bills).cwExp};
cwexp = cwexp(inds);

trenloc = cell(length(thesemice),numdates);
testloc = cell(length(thesemice),1);
for iMouse = 1:length(thesemice)
    if ~strcmp(datalocs(iMouse),'old') % dont bother with data on zserver
        
        subjectDates = {db(contains({db(:).mouse_name},thesemice{iMouse})).date};
        earliest_recording = subjectDates{1};
        
       
        [~,~,~,~,blkDir] = expDirs(thesemice{iMouse},earliest_recording,[],datalocs{iMouse});
        inds = strfind(blkDir,filesep);
        pardir = blkDir(1:inds(end-1));
        foo = dir(pardir);
        alldates = {foo.name};
        dis = find(contains(alldates,earliest_recording));
        whichDates = rectify(dis-numdates:dis-1,1);
        trenDates = alldates(whichDates);
        for iDate = 1:length(trenDates)
            [~, ~,hasBlock,pars] = dat.whichExpNums(thesemice{iMouse}, trenDates{iDate});
            pars = [pars{hasBlock}];
            x = [pars(:).distBetweenTargets];
            y = [pars(:).targetAltitude];
            trenloc{iMouse,iDate} = [x;y];
        end
        [~,~,~,pars] = dat.whichExpNums(thesemice{iMouse}, earliest_recording);
        p = pars{cwexp{iMouse}};
        x = [p(:).distBetweenTargets];
        y = [p(:).targetAltitude];
        testloc{iMouse} = [x;y];
    end
end

medianTrainingPos = nan(length(bills),2);
finalTrainingPos = nan(length(bills),2);
firstTestPos = nan(length(bills),2);
for iClu = 1:length(bills)
    k = kays(iClu);
    if ~strcmp(datalocs(k),'old') 
        medianTrainingPos(iClu,:) = median([trenloc{k,:}],2);
        finalTrainingPos(iClu,:) = median(trenloc{k,end},2);
        firstTestPos(iClu,:) = testloc{iMouse};
    end
end


%% PROFIT
[X,Y] = meshgrid(snrf.XPos,snrf.YPos);
allPos = [X(:) Y(:)];
mirrorPos = [-X(:) Y(:)];

normMain = nan(length(cfid),1);
normAnti = nan(length(cfid),1);
baba = nan(length(cfid),1);
for iClu = cfid'
    rf = clu_fits.RF{cfid == iClu};
    feet = clu_fits.bestFit{cfid == iClu}(:,:,1);
    
    thisPlace = mainCenter(cfid == iClu,:);
    dd_main = allPos - thisPlace;
    dd_anti = mirrorPos - thisPlace;
    mainLoc = sqrt(diag(dd_main*dd_main')) <= std_main(cfid == iClu);
    antiLoc = sqrt(diag(dd_anti*dd_anti')) <= std_main(cfid == iClu);
    normMain(cfid == iClu) = norm(rf(mainLoc),2);
    normAnti(cfid == iClu) = norm(rf(antiLoc),2);
end



