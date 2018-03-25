bs = 0.01;
tbin = 0.2;

rc = unique(right_cont);
lc = unique(left_cont);
resp = zeros(length(lc),length(rc),nCells);
errs = zeros(length(lc),length(rc),nCells);
psth = zeros(length(lc)*length(rc),round(tbin/bs + 1),nCells);
rast = cell(length(lc),length(rc),nCells);
for iclu = 1:nCells
    nt = zeros(length(lc),length(rc));
    for ll = 1:length(lc)
        for rr = 1:length(rc)
            theseStims = stimTimes(left_cont == lc(ll) & right_cont == rc(rr));
            nt(ll,rr) = length(theseStims);
            [p,bins,rx,ry,nspk] =  ... 
                psthRasterAndCounts(spks.st(spks.clu == inclCID(iclu)),theseStims,[0 tbin],bs);
            [~,~,~,~,base] =  ... 
                psthRasterAndCounts(spks.st(spks.clu == inclCID(iclu)),theseStims,[-0.05 0.05],0.01);
            psth(ll + 4*(rr-1),:,iclu) = p;
            resp(ll,rr,iclu) = mean(nspk-base);
            errs(ll,rr,iclu) = std(nspk-base);
            rast{ll,rr,iclu} = [rx;ry];
        end
    end
end
nr = sum(nt,1)/max(sum(nt,1))+0.5;

%%
cwstart = aln.cw2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
cwend = max(stimTimes(iscw));
passtart = aln.pas2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
pasend = max(stimTimes(~iscw));
figure('units','normalized','position',[0.0370 0.2231 0.8786 0.6389]);
for iclu = 1:nCells
    dis = find(spks.cids == inclCID(iclu));
    subplot(2,4,1) % psth
    plot(bins,mean(psth(1:4,:,iclu),1),'linewidth',nr(1))
    hold on;
    plot(bins,mean(psth(5:8,:,iclu),1),'linewidth',nr(2))
    plot(bins,mean(psth(9:12,:,iclu),1),'linewidth',nr(3))
    plot(bins,mean(psth(13:16,:,iclu),1),'linewidth',nr(4))
    legend(split(num2str(lc')))
    ylabel('Hz');
    hold off;
    subplot(2,4,2) % contrast tuning curve
    plot(rc,mean(resp(:,:,iclu),1))
    ylabel('mean response')
    hold off
    subplot(2,4,5) % all psth
    imagesc(bins,[],psth(:,:,iclu))
    hold on;
    plot([xlim NaN xlim NaN xlim],[4.5 4.5 NaN 8.5 8.5 NaN 12.5 12.5],'--','color','w','linewidth',2)
    tix = yticks;
    lab = repmat(lc,4,1);
    set(gca,'YTickLabels',lab(tix))
    ylabel('left cont')
    xlabel('time')
    hold off;
    subplot(2,4,6) % full contrast tuning
    imagesc_but_good(rc,lc,resp(:,:,iclu))
    hold on;
    text(repelem(1:4,1,4),repmat(1:4,1,4),split(num2str(reshape(nt,16,1)')))
    xlabel('right'); ylabel('left')
    hold off;
    
    subplot(2,2,2) % amplitude 
    scatter(spks.st(spks.clu == inclCID(iclu)),spks.gain*spks.spikeAmps(spks.clu == inclCID(iclu)),'.')
    hold on;
    linesAtEvents([cwstart cwend],gca,'--','linewidth',2,'color','b')
    linesAtEvents([passtart pasend],gca,'--','linewidth',2,'color','r')
    ylabel('Amplitude')
    title(['Cell ' num2str(inclCID(iclu)) ', #spks = ' num2str(sum(spks.clu == inclCID(iclu)))])
    hold off
    subplot(2,2,4) % receptive field
    imagesc(snrf.XPos,snrf.YPos,snrf.neur_rfmap{dis})
    set(gca,'ydir','normal')
    hold on; 
    scatter(snrf.neur_rfstats(dis).fit2D(2),snrf.neur_rfstats(dis).fit2D(3),120,'k','*')
    plot(cwSigma*cos(0:0.1:2*pi) + cwXPos,cwSigma*sin(0:0.1:2*pi) + cwYPos,'--','color','w','linewidth',2)
    [~, ind1] = max(max(abs(snrf.neur_rfmap{dis}),[],2));
    [~, ind2] = max(max(abs(snrf.neur_rfmap{dis}),[],1));
    moddir = sign(snrf.neur_rfmap{dis}(ind1,ind2));
    title(moddir*snrf.neur_rfstats(dis).peakZscore)
    hold off;
    
    pause
end
