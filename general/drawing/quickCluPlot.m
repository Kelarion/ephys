function quickCluPlot(spks,CID)
% quickCluPlot(spks,CID)
%

thresh = 0.2;

gcf;
subtightplot(4,3,[1 2 3],[],[0.01 0.01],[0.1 0.01])
scatter(spks.st(spks.clu == CID),spks.gain*spks.spikeAmps(spks.clu == CID),24,'markeredgealpha',0.3)
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

subtightplot(4,3,[5 6 8 9 11 12],[0.5,0.01]); % waveform
clut = spks.cluTemps(CID+1)+1;
[~,ind] = max(max(spks.tempsUnW(clut,:,:)));
wf = permute(spks.tempsUnW(clut,:,:),[3 2 1]);

chanAmps = max(wf,[],2)-min(wf, [], 2);
maxAmp = max(chanAmps);
inclchans = chanAmps>maxAmp*thresh;

plotchans = rectify(ind-9,1):saturate(ind+9,374);

plotWaveform(wf(plotchans,:),spks.xcoords(plotchans),spks.ycoords(plotchans),20,5,0,'k')
hold on
plotWaveform(wf(inclchans,:),spks.xcoords(inclchans),spks.ycoords(inclchans),20,5,0,'b')
hold off
ylabel('Depth from bottom')
title(sprintf('Cell %d template and amplitude',CID))


end