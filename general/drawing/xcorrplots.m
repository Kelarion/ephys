nchan = 10;
istart = 11;

goodch = spks.cids(spks.cgs == 2);
chans = goodch(istart:(istart+nchan));
j = 1;
for row = 1:nchan
    for col = 1:nchan
        if row > col
            j = j+1;
            continue
        end
        ax = subplot(nchan,nchan,j);
        hold off;
        if row == col
            myACG(spks.st(spks.clu == chans(row)), ax, []);
            ylabel(['cluster ' num2str(chans(row))])
            ax.XTick = [];
        else
            spikesXCG(spks.st(spks.clu == chans(row)), spks.st(spks.clu == chans(col)),ax);
            ax.XTick = [];  ax.YTick = [];
        end
        if row == 1
            title(['cluster ' num2str(chans(col))])
        end
        j = j+1;
    end
end
