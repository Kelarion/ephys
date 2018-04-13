binsize = 0.001; % bin size
win = [0 0.2]; % maximum time

rc = unique(right_cont);
lc = unique(left_cont);
resp = zeros(length(lc),length(rc),nCells);
errs = zeros(length(lc),length(rc),nCells);
psth = zeros(length(lc)*length(rc),round(diff(win)/binsize),nCells);
for iclu = 1:nCells
    nt = zeros(length(lc),length(rc));
    [~,bins,~,~,~,ba] = psthAndBA(spks.st(spks.clu == inclCID(iclu)),stimTimes,win,binsize);
    for ll = 1:length(lc)
        for rr = 1:length(rc)
            theseStims = left_cont == lc(ll) & right_cont == rc(rr);
            p = mean(ba(theseStims,:)./binsize,1);
            
            nt(ll,rr) = sum(theseStims);
            resp(ll,rr,iclu) = mean(sum(ba(theseStims,:),2)/0.2);
            errs(ll,rr,iclu) = std(sum(ba(theseStims,:),2)/0.2)./sqrt(nt(ll,rr));
            psth(ll + 4*(rr-1),:,iclu) = p;
        end
    end
end
nr = sum(nt,1)/max(sum(nt,1))+0.5;

%%
cwstart = aln.cw2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
cwend = max(stimTimes(iscw));
passtart = aln.pas2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
pasend = max(stimTimes(~iscw));
% figure('units','normalized','position',[0.0370 0.2231 0.8786 0.6389]);
for iclu = 1:nCells
    [resp,errs,nt,whichEvs,psth,bins,rast] = tuningCurve2D(spks.st(spks.clu == inclCID(iclu)), ...
        stimTimes,left_cont,right_cont,win,binsize);
    
    plotTuningCurveTemp
    
end
