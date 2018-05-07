%% ACGs
nw = 1.5;
figure;
for ii = 1:length(feats.nSpks)
    t = feats.ACG_bins(1,:);
    ACG = feats.ACG(ii,:);
    smoth = smooth(ACG,0.02,'loess');

    subplot(1,3,1)
    plot(t,ACG);
    hold on; 
    plot(t,smooth(ACG,10),'linewidth',3);
    t_m = dot(t(1:200),ACG(1:200))/sum(ACG(1:200));
    r_m = mean(ACG);
    plot([t_m t_m],ylim,'--','color','k');
    hold off;
    title(['Cell: ' num2str(feats.whichCell(ii,1)) ' dset: ' num2str(feats.whichCell(ii,3))])
    ylabel('Conditional intensity (Hz)')
    xlabel('Lag (sec)')
    
    subplot(1,3,2)
    plot(feats.F(1:10:end),feats.powerSpectra(ii,1:10:end))
    [pmax, ind] = max(feats.powerSpectra(ii,:));
    hold on
    plot(feats.F,smooth(feats.powerSpectra(ii,:),0.01,'loess'),'linewidth',3)
    plot(feats.F([ind ind]),ylim,'--')
    hold off
    xlabel('frequency (Hz)')
    ylabel('power ? (dB???)')
%     tsp = (1:size(feats.spikeShapes,2))/30000;
%     spwf = (feats.spikeShapes(ii,:)/range(feats.spikeShapes(ii,:)))*feats.amplitude(ii);
%     plot(tsp,spwf,'linewidth',2)
    title([num2str(pmax) ' at ' num2str(feats.F(ind)) 'Hz'])
    
    subplot(1,3,3)
    plot(feats.spikeShapes(ii,:))
    title(num2str(feats.amplitude(ii)))
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

%% return map (joint interval distributions)
figure; 
for ii = 1:length(spks.cids)
    myisi = diff(spks.st(spks.clu == spks.cids(ii)));
    nspk = sum(spks.clu == spks.cids(ii));
    pad = nan(abs(mod(nspk,2)-1),1);
    [previsi, ord] = sort(myisi(1:2:end));
    previsi = previsi*1000;
    subsisi = [myisi(2:2:end); pad];
    subsisi = subsisi(ord)*1000;
    subplot(1,2,1)
    scatter(previsi,subsisi,'.','markeredgecolor',[0.4 0.4 0.4],'markeredgealpha',0.3)
    hold on; set(gca,'xscale','log');set(gca,'yscale','log')
    title(['Cell: ' num2str(spks.cids(ii))])
    ylabel('Next ISI, i.e. 2nd-order')
    xlabel('ISI (ms)')
    hold off
    shiftedIsi = [];
    for tau = 1:40; shiftedIsi = [shiftedIsi myisi(tau:40:end-(41-tau))];end
    sercor = corrcoef(shiftedIsi);
    subplot(1,2,2)
    plot(2:40,sercor(1,2:end))
    set(gca,'ylim',[-1 1])
    xlabel('Order')
    ylabel('Correlation coef')
    pause
end

%% spikes per event
nbin = length(feats.ISI_bins);
thr = 0.006;
figure
for iCell = 7:nClu
    thesest = spks.st(spks.clu == spks.clu(iCell));
    isi = diff(thesest);
    isiBin = feats.ISI_bins;
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


