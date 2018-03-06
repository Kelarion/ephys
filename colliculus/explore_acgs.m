%%
figure;
for ii = 1:967
    t = feats.ACG_bins(ii,:);
    ACG = feats.ACG(ii,:);
    smoth = smooth(ACG,0.02,'loess');
    
    subplot(1,2,1)
    plot(t,ACG);
    hold on; 
    plot(t,smooth(ACG,10),'linewidth',3);
    t_m = dot(t(1:200),ACG(1:200))/sum(ACG(1:200));
    r_m = mean(acg);
    plot([t_m t_m],ylim,'--','color','k');
    hold off;
    title(['Num spikes: ' num2str(feats.nSpks(ii))]);
    ylabel('Conditional intensity')
    xlabel('Lag (sec)')
    
    subplot(1,2,2)
    tsp = (1:size(feats.spikeShapes,2))/30000;
    spwf = (feats.spikeShapes(ii,:)/range(feats.spikeShapes(ii,:)))*feats.amplitude(ii);
    plot(tsp,spwf,'linewidth',2)
    title(['Cell: ' num2str(feats.whichCell(ii,1)) ' dset: ' num2str(feats.whichCell(ii,2))])
    pause;
end

%%
figure;
toy = 'a*exp(-b*x)';
for ii = 1:967
    [m,ind] = max(feats.distISI(ii,:));
    t = feats.ISI_bins(ind:ind+200);
    pisi = feats.distISI(ii,ind:ind+200)/sum(feats.distISI(ii,ind:ind+200));
    smoth = smooth(pisi,0.02,'loess');
    subplot(1,2,1)
    scatter(t,pisi);
    set(gca,'yscale','log')
    hold on; 
%     plot(t,smooth(pisi,10),'linewidth',3);
    t_m = dot(t,pisi)/sum(pisi);
    t_var = dot((t-t_m).^2,pisi)/sum(pisi);
    f = fit(t',pisi',toy,'StartPoint',[m 1/t_m]);
    plot(t,f(t),'linewidth',2,'color','r')
    plot([t_m t_m],ylim,'--','color','k');
    hold off;
    title(['Num spikes: ' num2str(feats.nSpks(ii))]);
    ylabel('Total count')
    xlabel('ISI (sec)')
    
    subplot(1,2,2)
    tsp = (1:size(feats.spikeShapes,2))/30000;
    spwf = (feats.spikeShapes(ii,:)/range(feats.spikeShapes(ii,:)))*feats.amplitude(ii);
    plot(tsp,spwf,'linewidth',2)
    title(['Cell: ' num2str(feats.whichCell(ii,1)) ' dset: ' num2str(feats.whichCell(ii,2))])
    
    pause;
end

%%
figure; 
for ii = 1:length(spks.cids)
    myisi = diff(spks.st(spks.clu == spks.cids(ii)));
    nspk = sum(spks.clu == spks.cids(ii));
    pad = nan(abs(mod(nspk,2)-1),1);
    [previsi, ord] = sort(myisi(1:2:end));
    previsi = previsi*1000;
    subsisi = [myisi(2:2:end); pad];
    subsisi = subsisi(ord)*1000;
    scatter(previsi,subsisi,'.','markeredgecolor',[0.4 0.4 0.4])
    hold on; set(gca,'xscale','log');set(gca,'yscale','log')
    title(num2str(spks.cids(ii)))
    hold off
    pause
end