function fig = behTraj(trajes, colpart, tLims)


t = [1:size(trajes,2)]./200;

if nargin < 3 || isempty(tLims)
    tLims = [t(1) t(end)];
end

fig = figure; 
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);

set(ax1,'Position',[0.1300 0.4358 0.7750 0.2157]);
plot(ax1,t,trajes(1,:),'Color',colpart,'LineWidth',2);

set(ax2,'Position',[0.1300 0.2869 0.7750 0.2157]);
plot(ax2,t,trajes(2,:),'Color',colpart,'LineWidth',2);

set(ax3,'Position',[0.1300 0.2048 0.7750 0.2157]);
plot(ax3,t,trajes(3,:),'Color',colpart,'LineWidth',2);

linkaxes(fig.Children); 
set(ax1,'XLim',tLims)
% axis(ax1,'off');axis(ax2,'off');axis(ax3,'off');

end