%% ACGS
nw = 1.5;
figure;
for ii = 1:length(feats.nSpks)
    t = feats.ACG_bins;
    ACG = feats.ACG(ii,:);
    smoth = smooth(ACG,0.02,'loess');
    
    subplot(1,2,1)
    plot(t,ACG);
    hold on; 
    plot(t,smooth(ACG,10),'linewidth',3);
    t_m = dot(t(1:200),ACG(1:200))/sum(ACG(1:200));
    r_m = mean(ACG);
    plot([t_m t_m],ylim,'--','color','k');
    hold off;
    title(['Cell: ' num2str(feats.whichCell(ii,1)) ' dset: ' num2str(feats.whichCell(ii,2))])
    ylabel('Conditional intensity')
    xlabel('Lag (sec)')
    
    subplot(1,2,2)
    [pxx,f] = pmtm(smooth(ACG,5) - mean(ACG(end-50:end)),nw,0:0.01:200,1/mean(diff(feats.ACG_bins)));
    [pmax, ind] = max(10*log10(pxx));
    plot(f,10*log10(pxx))
%     set(gca,'yscale','log')
    hold on; plot(f([ind ind]),ylim,'--')
    hold off
    xlabel('frequency (Hz)')
%     tsp = (1:size(feats.spikeShapes,2))/30000;
%     spwf = (feats.spikeShapes(ii,:)/range(feats.spikeShapes(ii,:)))*feats.amplitude(ii);
%     plot(tsp,spwf,'linewidth',2)
    title([num2str(pmax) 'dB at ' num2str(f(ind)) 'Hz'])
    pause;
end

%% isi distribution with exponential fit
figure;
toy = 'a*x + b';
for ii = 1:size(feats.features,1)
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
    f = fit(t',log10(pisi+0.00001)',toy,'StartPoint',[diff(pisi([1 end])) log10(m)]);
    plot(t,10.^(f(t)),'linewidth',2,'color','r')
%     plot([t_m t_m],ylim,'--','color','k');
    plot([0.06 0.06],ylim,'--','color','k')
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

%% return map
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

%% spikes per event
nbin = isi_time/isi_binsize;
thr = 0.006;
figure
for iCell = 7:nClu
    thesest = spks.st(spks.clu == goodMesoCIDs(iCell));
    isi = diff(thesest);
    isiBin = 0:isi_binsize:isi_time;
    p = histcounts(isi,isiBin);
    
    burstInds = find(isi<thr);
    ibi = [diff(burstInds); 2];
    nEv = zeros(sum(isi<thr),1);
    ii = 1;
    while ii < sum(isi<thr)
        nEv(ii) = find(ibi(ii:end)>1,1,'first');
        ii = ii + nEv(ii);
    end
    burstInds = burstInds(nEv > 0);
    spkPerBurst = nEv(nEv>0);
    spb = 1:max(spkPerBurst);
    if sum(nEv>0) < 1, continue; end
    pNspk = histcounts(spkPerBurst,length(spb));
    pNspk = pNspk/sum(pNspk);
    
    egBurst = burstInds(find(spkPerBurst == max(spkPerBurst),1,'first'));
    
    subplot(2,2,3)
    stairs(isiBin(isiBin<=0.2),p(isiBin<=0.2))
    hold on; plot([thr thr],ylim,'--','color','k')
    hold off
    subplot(2,2,4)
    scatter(spb,log10(pNspk))
    ugh = log10(pNspk(:)) > -Inf;
    try
        f = fit(spb(ugh)',log10(pNspk(ugh))','a*x + b','startpoint',[-1 0]);
        hold on; plot(spb,f(spb),'--','color','k');
    catch
    end
    hold off
    
    subplot(2,1,1)
    rasterdot(thesest((egBurst-6):(egBurst+max(spkPerBurst)+6)));
    title(num2str(goodMesoCIDs(iCell)))
    pause
end

%%
figure;
winsize = 5; % mins
nwin = round(max(thesest)/(60*winsize));
j = 1;
for ii = 1:nwin
    fin = find(thesest < thesest(j)+60*winsize,1,'last');
    tic
    [S, f, R] = mtspectrumpt(thesest(j:fin),params);
    plot(f,S)
    toc
    pause
end


