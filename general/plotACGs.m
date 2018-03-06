function plotACGs(ACGs,ACG_bins,iStart,nchan)
% plotACGs(ACGs,ACG_bins,iStart,nchan)



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
            ACGs;
            ylabel(['cluster ' num2str(chans(row))])
            ax.XTick = [];
        else
            stairs(lags,n,'linewidth',1.5)
            hold on; plot([0 0],ax.YLim,'--','color','k'); hold off
            ax.XTick = [];  ax.YTick = [];
        end
        if row == 1
            title(['cluster ' num2str(chans(col))])
        end
        j = j+1;
    end
end

end