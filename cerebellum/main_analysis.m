%% load 
load('C:\Users\Matteo A\Documents\MATLAB\mine\ephys\cerebellum\allrips.mat')
ksRoot = 'I:\data\';
localKsRoot = 'C:\DATA\Spikes\Janelia\';
saveDir = 'C:\DATA\Spikes\Janelia\Videos\';

% To convert from the 3-number identification (which correspond to cell ID,
% tag, and dataset) to a single unique cell ID, I use the below mappings. 
% They're similar to a conversion between base 5 and 10. For some 
% reason, using base-5 makes this mapping always invertible. (I checked 
% this with randomly-generated ID matrices where the 'digits' could vary 
% between 1 and 2000.)
b = 5;
wc2id = @(WC) WC*[b^2 b^1 b^0]';
id2wc = @(id) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)];

CID = wc2id(R.whichCell);
cgs = unique(CID);
nClu = length(cgs);
%% assess significance of phase-locking 
% in reality, testing for non-uniformity of angle distribution
minNSpk = 50;
alph = 0.05;

pvals = nan(nClu,3); 
for iClu = 1:nClu
    if (sum(CID == cgs(iClu)) > minNSpk)
        pvals(iClu,1) = circ_rtest(R.spkPhase(CID == cgs(iClu)));
    end
    if (sum(CID == cgs(iClu) & ~R.isOpto) > minNSpk)
        pvals(iClu,2) = circ_rtest(R.spkPhase(CID == cgs(iClu) & ~R.isOpto));
    end
    if (sum(CID == cgs(iClu) & R.isOpto) > minNSpk)
        pvals(iClu,3) = circ_rtest(R.spkPhase(CID == cgs(iClu) & R.isOpto));
    end
end

islocked = sum(pvals(:,[1 2]) < alph,2) == 2; % we care about phase-locking spontaneously
islocked = islocked & (pvals(:,3) < alph | isnan(pvals(:,3)));

%% mean phase comparisons
[cidOp,~,indOp] = unique(CID(R.isOpto));
[optoPhs,optoMag] = splitapply(@circmean,R.spkPhase(R.isOpto),indOp);

[cidSp,~,indSp] = unique(CID(~R.isOpto));
[spontPhs,spontMag] = splitapply(@circmean,R.spkPhase(~R.isOpto),indSp);
spontFreq = splitapply(@mean,R.freq(~R.isOpto),indSp);
spontPk = splitapply(@mean,R.magnitude(~R.isOpto),indSp);

inboth = intersect(cidSp,cidOp); % which cells fired in both events
spNop = ismember(cidSp,inboth);
opNsp = ismember(cidOp,inboth); % get indices of these cells
spontPhs = spontPhs(spNop);
spontMag = spontMag(spNop);
spontFreq = spontFreq(spNop);
spontPk = spontPk(spNop);
optoPhs = optoPhs(opNsp);
optoMag = optoMag(opNsp);

lckd = ismember(inboth,cgs(islocked));

mags = spontMag; % determine the magnitude of each point

dtheta = spontPhs - optoPhs; % plot 
phsDiff = mags.*exp(1i*dtheta);
phsors = reshape([zeros(size(phsDiff)) phsDiff nan(size(phsDiff))].',length(inboth)*3,1);

% plot
figure
subplot(1,2,1)
scatlab(optoPhs(lckd),spontPhs(lckd),inboth(lckd),64*mags(lckd),spontFreq(lckd),'filled')
hold on
scatlab(optoPhs(~lckd),spontPhs(~lckd),inboth(~lckd),48*mags(~lckd),spontFreq(~lckd))
samelims
xlabel('phase during induced')
ylabel('phase during spontaneous')
title('Mean phase-locking')
grid on; box on
plot(xlim,xlim,'--','color','k')
plot(xlim,xlim - 2*pi,'--','color','k')
plot(xlim,xlim + 2*pi,'--','color','k')
hold off

subplot(1,2,2)
set(gca,'XLim',[-1, 1])
set(gca,'YLim',[-1, 1])
hold on
% make circular grid
for rr = [0.25 0.5 0.75 1]
    plot(rr*sin(0:0.1:2.1*pi),rr*cos(0:0.1:2.1*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
for rr = [0.25 0.5 0.75 1]'
    plot(rr*sin(0:pi/8:2*pi),rr*cos(0:pi/8:2*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
% ------------------
plot(real(phsors),imag(phsors),'linewidth',2)
scatlab(real(phsDiff(lckd)),imag(phsDiff(lckd)),inboth(lckd),[],spontFreq(lckd),'filled')
scatlab(real(phsDiff(~lckd)),imag(phsDiff(~lckd)),inboth(~lckd),[],spontFreq(~lckd))

title('\boldmath${\bar{z}_{spont} / \hat{\bar{z}}_{opto}}$','fontsize',14,'interpreter','latex')
xlabel('Real part');ylabel('Imaginary part')

c = colorbar;
c.Label.String = 'mean spontaneous frequency (Hz)';
c.Label.FontSize = 11;
% colormap jet

%% -----------------------------
foo = cgs(islocked);
eg_cid = foo(randi(length(foo)));

wrapped_spont = exp(1i*(R.spkPhase(CID == eg_cid & ~R.isOpto)));
[theta, r] = circmean(R.spkPhase(CID == eg_cid & ~R.isOpto));
z_spont = r*exp(1i*theta);

wrapped_opto = exp(1i*(R.spkPhase(CID == eg_cid & R.isOpto)));
[theta, r] = circmean(R.spkPhase(CID == eg_cid & R.isOpto));
z_opto = r*exp(1i*theta);

figure
subplot(1,2,1);
xlim([-1 1]);ylim([-1 1]);
for rr = [0.25 0.5 0.75 1]
    plot(rr*sin(0:0.1:2.1*pi),rr*cos(0:0.1:2.1*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
    hold on
end
for rr = [0.25 0.5 0.75 1]'
    plot(rr*sin(0:pi/8:2*pi),rr*cos(0:pi/8:2*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
scatter(real(wrapped_spont),imag(wrapped_spont),'.')
plot([0 real(z_spont)],[0 imag(z_spont)],'linewidth',2)
scatter(real(z_spont),imag(z_spont),'r','filled')
% axis off
hold off

subplot(1,2,2);
xlim([-1 1]);ylim([-1 1]);
for rr = [0.25 0.5 0.75 1]
plot(rr*sin(0:0.1:2.1*pi),rr*cos(0:0.1:2.1*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
hold on
end
for rr = [0.25 0.5 0.75 1]'
plot(rr*sin(0:pi/8:2*pi),rr*cos(0:pi/8:2*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
scatter(real(wrapped_opto),imag(wrapped_opto),'.')
plot([0 real(z_opto)],[0 imag(z_opto)],'linewidth',2)
scatter(real(z_opto),imag(z_opto),'r','filled')
% axis off
hold off


%% timing of events
interval_bins = linspace(0,10,51); % inter-event intervals
oid = unique(wc2id(R.osc.whichCell));

start2start = zeros(length(interval_bins)-1,length(oid));
end2end = zeros(length(interval_bins)-1,length(oid));
end2start = zeros(length(interval_bins)-1,length(oid));

for iEv = 1:length(oid)
    evTimes = R.osc.Times(wc2id(R.osc.whichCell) == oid(iEv) & isnan(R.osc.optoFreq),:);
    start2start(:,iEv) = histcounts(diff(evTimes(:,1)),interval_bins);
    end2end(:,iEv) = histcounts(diff(evTimes(:,2)),interval_bins);
    end2start(:,iEv) = histcounts(diff([evTimes(2:end,1) evTimes(1:end-1,2)]),interval_bins);
end

% % for only looking at the phase-locked events
% start2start = zeros(length(interval_bins)-1,sum(islocked));
% end2end = zeros(length(interval_bins)-1,sum(islocked));
% end2start = zeros(length(interval_bins)-1,sum(islocked));
% 
% cids = cgs(islocked);
% for iClu = 1:sum(islocked)
%     cluTimes = R.interval(CID == cids(iClu) & ~R.isOpto,:);
%     [~,theseEvs] = unique(R.whichEvent(CID == cids(iClu) & ~R.isOpto));
%     start2start(:,iClu) = histcounts(diff(cluTimes(theseEvs,1)),interval_bins);
%     end2end(:,iClu) = histcounts(diff(cluTimes(theseEvs,2)),interval_bins);
%     foo = cluTimes(theseEvs,:);
%     end2start(:,iClu) = histcounts(diff([foo(2:end,1) foo(1:end-1,2)]),interval_bins);
% end

figure;
bar(interval_bins(1:end-1),sum(end2start,2),'histc')
xlabel('Inter-event interval (sec)')
ylabel('Number of events')

% --------------------

durbins = linspace(4,45,41); % spontaneous event duration distribution
n = histcounts(diff(R.osc.Times(isnan(R.osc.optoFreq),:),[],2).*R.osc.Freq(isnan(R.osc.optoFreq)),durbins,...
    'normalization','probability');

durbins = durbins(1:end-1);
linmod = polyfit(durbins(n>0),log10(n(n>0)),1);

figure;
subplot(2,1,1)
bar(durbins,n,'histc')
ylabel('Probability')
xticks([])

subplot(2,1,2)
scatter(durbins,n,'filled')
hold on
plot(xlim,10.^([xlim; 1 1]'*linmod'))
set(gca,'yscale','log')
ylabel('log10(Probability)')
xlabel('Number of cycles in event')
text(max(xlim),max(ylim),['p(cycle) = ' num2str(10^linmod(1))], ... 
    'horizontalalignment','right','verticalalignment','top')

