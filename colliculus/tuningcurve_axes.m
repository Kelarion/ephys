% Makes the axis layout for plotting the 2D tuning curves, as well as some
% additional useful information for each cluster

% Needs: binned array of spikes, spks, snrf, cwSigma, cwXPos, cwYPos,
% cwstart, cwend, passtart, pasend, inclCID, iclu


%% tuning surface with marginals
ax(1) = subtightplot(5,9,[10 11 19 20]); % full tuning surface

ax(2) = subtightplot(5,9,[1 2]); % upper marginal

ax(3) = subtightplot(5,9,[12 21]); % right marginal

%% receptive field
ax(4) =subtightplot(5,9,[28 29 30 37 38 39],[0.15 0.01]);
imagesc(snrf.XPos,snrf.YPos,snrf.neur_rfmap{snrf.neur_ID == inclCID(iclu)})
set(gca,'ydir','normal')
hold on;
plot(1.5*cwSigma*cos(-1:0.1:2*pi) + cwXPos,1.5*cwSigma*sin(-1:0.1:2*pi) + cwYPos,'-k','linewidth',2)
plot(1.5*cwSigma*cos(-1:0.1:2*pi) - cwXPos,1.5*cwSigma*sin(-1:0.1:2*pi) + cwYPos,'-k','linewidth',2)
title('Sparse noise RF')
hold off;

%% psth
ax(5) =subtightplot(5,9,[ 4 5 6],[0.01 0.1],[0.01, 0.05],[0.05 0.05]); % mean psth

ax(6) =subtightplot(5,9,[ 13 14 15  22 23 24],[0.01 0.1]); % all psth

ax(7) =subtightplot(5,9,[ 31 32 33  40 41 42],[0.01 0.1]); % raster
[tr,b] = find(ba);
[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3);
rasterScale = 6;
rasterY(2:3:end) = rasterY(2:3:end)+rasterScale;
plot(rasterX,rasterY, 'k','linewidth',2);
ylim([min(rasterY) max(rasterY)])
hold(ax(7),'on')
plot(ax(7),[0 0],ylim(ax(7)),'-','color','r','linewidth',2)
box off
hold(ax(7),'off')
linkaxes(ax(5:7),'x')
ylabel('Event number')

%% cluster information
ax(8) =subtightplot(5,9,[7 8 9],[0.05,0.05]); % amplitude plot
scatter(spks.st(spks.clu == inclCID(iclu)),spks.gain*spks.spikeAmps(spks.clu == inclCID(iclu)),24,'markeredgealpha',0.1)
% hold on
% linesAtEvents([cwstart cwend],gca,'--','linewidth',2,'color','b')
% text(mean([cwstart cwend]),max(ylim),'active','color','b','fontweight','bold','fontsize',14, ...
%     'verticalalignment','top','horizontalalignment','center')
% linesAtEvents([passtart pasend],gca,'--','linewidth',2,'color','r')
% text(mean([passtart pasend]),max(ylim),'passive','color','r','fontweight','bold','fontsize',14,...
%     'verticalalignment','top','horizontalalignment','center')
% hold off
ylabel('spike amp')
xlabel('time')

ax(9) =subtightplot(5,9,[16 17 43 44],[0.2,0.05]); % waveform
clut = spks.cluTemps(inclCID(iclu)+1)+1;
[~,ind] = max(max(spks.tempsUnW(clut,:,:)));
plotchans = rectify(ind-7,1):saturate(ind+7,374);
wf = permute(spks.tempsUnW(clut,:,plotchans),[3 2 1]);
hold off
plotWaveform(wf,spks.xcoords(plotchans),spks.ycoords(plotchans),20,5,0,'k')
ylabel('Depth from bottom')
title('Spike template')

ax(10) =subtightplot(5,9,[18 27 36 45],[0.1,0.01]); % depth
xlim([-0.5 0.5])
plot(xlim,[scTop scTop],'k','linewidth',2)
text(-0.5,scTop,' SC surface','verticalalignment','bottom')
hold on;
plot(xlim,[scBottom scBottom],'k','linewidth',2)
text(-0.5,scBottom,' sSC/dSC','verticalalignment','top')
ylim([0 3840])
text(0,mean(spks.spikeDepths(spks.clu == inclCID(iclu))),'X', ... 
    'fontweight','bold','fontname','arial','fontsize',14,...
    'verticalalignment','middle');
xticks([])
set(gca,'yaxislocation','right')
hold off;

