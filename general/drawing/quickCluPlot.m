function quickCluPlot(spks,CID)
% quickCluPlot(spks,CID)
%

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
plotchans = rectify(ind-7,1):saturate(ind+7,374);
wf = permute(spks.tempsUnW(clut,:,plotchans),[3 2 1]);
plotWaveform(wf,spks.xcoords(plotchans),spks.ycoords(plotchans),20,5,0,'k')
ylabel('Depth from bottom')
title('Spike template')


end