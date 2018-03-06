function [n, lags] = spikesXCG(st1, st2, ax)
% [n, lags] = spikesXCG(st1, st2, ax)
% CCG centred on st1

binSize = 0.01;

b = -0.5:binSize:0.5;

[n, lags] = histdiff(st2, st1, b);
n = n./binSize; 
n = n(:); lags = lags(:); 
n = n/(binSize*length([st1; st2])); % divide by P(spike2) and convert to Hz

if ~isempty(ax)
    stairs(lags,n,'linewidth',1.5)
    hold on; plot([0 0],ax.YLim,'--','color','k'); hold off
end

