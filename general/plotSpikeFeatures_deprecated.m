function Ax = plotSpikeFeatures(feets,ptype,Fs,whichAP)
% plotSpikeFeatures(feets,ptype,Fs,whichAP)
%
% makes scatterplot matrix of the spike waveform features calculated by
% the 'spikeFeatures' function, which makes the first input structure
%
% Inputs:
%   - feets [struct]: output of spikeFeatures
%   - ptype: 'scatter' to make regular scatterplots 
%            'density' to plot the density (default)
%   - Fs: sampling rate (default to imec: 30kHz)
%   - whichAP (optional): which spike to use as illustrative example, will
%            pick a random one by default
%
% to do: make interactive, update example spike by clicking a region of the
% scatter plot

if nargin<4, whichAP = randperm(size(feets.whichCell,1),1); end
if nargin<3, Fs = 30000; end
if nargin<2, ptype = 'density'; end

egAP = smooth(feets.spikeShapes(whichAP,:),0.1,'loess');
[t, tTro] = min(egAP);
[p, tPeak] = max(egAP(tTro:end));
tPeak = tPeak + tTro - 1;

% re-compute things for plotting
k = 1/(Fs/1000);
t_temp = [1:length(egAP)]*k;
hwhm = (feets.FWXM(whichAP,1)/2)*Fs;
tHM1 = (tTro - hwhm)*k; tHM2 = (tTro + hwhm)*k;

hw3m = (feets.FWXM(whichAP,2)/2)*Fs;
t3M1 = (tTro - hw3m)*k; t3M2 = (tTro + hw3m)*k;

G1 = feets.Grads(whichAP,1)/k; G2 = feets.Grads(whichAP,2)/k;
tG1 = linspace(tTro+0.8,tTro+3.2,4)*k; 
tG2 = linspace(tTro+9.5,tTro+14.5,4)*k;
tan1 = (tG1 - t_temp(tTro+2))*G1 + egAP(tTro+2); % drawing tangent lines
tan2 = (tG2 - t_temp(tTro+12))*G2 + egAP(tTro+12);

% main plot
allFeet = [feets.TP*1000 feets.FWXM*1000 feets.Grads feets.HRatio];
whichFeet = {'TtP','FWHM','FW3M','dV(0.06)','dV(0.4)','Peak/Trough'};
% niceWFs = isnan(feets.I_cap); % ignore cells with prominent capacitative phase
niceWFs = ~isnan(feets.TP); % ignore cells without prominent K+ current
nVar = size(allFeet,2);
j = 1;
for row = 1:nVar 
    for col = 1:nVar
        if row < col
            j = j+1;
            continue
        elseif row == col
            Ax(row,col) = subplot(nVar,nVar,j);
            histogram(allFeet(:,row),25,'edgealpha',0);
            Ax(row,col).YTick = [];
            box(Ax(row,col),'off')
            if col == 1, ylabel(whichFeet{row},'fontweight','bold'); end
            if row == nVar, xlabel(whichFeet{col},'fontweight','bold'); end
            j = j+1; 
            continue
        end
        Ax(row,col) = subplot(nVar,nVar,j);
        
        if strcmp(ptype,'density')
            [n, c] = hist3(allFeet(:,[col, row]),[15 15]);
            imagesc(c{1},flipud(c{2}(:)),flipud(n')); colormap gray;
            set(Ax(row,col),'Ydir','normal')
        else
            scatter(allFeet(niceWFs,col),allFeet(niceWFs,row),'.')
        end
        Ax(row,col).XTick = [0 Ax(row,col).XLim(end)]; 
        Ax(row,col).YTick = [0 Ax(row,col).YLim(end)];
        box(Ax(row,col),'on')
        
        if col == 1, ylabel(whichFeet{row},'fontweight','bold'); end
        if row == nVar, xlabel(whichFeet{col},'fontweight','bold'); end
        
        j = j+1;
    end
end
for ii = 1:nVar, linkaxes(Ax(:,ii),'x'); end

% add example waveform 
ax = subplot(2,2,2);
pos = ax.Position;

plot(t_temp,egAP,'linewidth',2)
ax.XLim = [tTro/2*k t_temp(end)];
% annotations are annoyingly in 'normalized' coordinates, so I need to
% convert from actual coordinates to values normalized for the figure, i.e.
% (0,0) is lower left, (1,1) upper right. 'pos' tells you the position of
% the axes in these 'normalized' coordinates.
x2norm = @(x)(x-ax.XLim(1))/diff(ax.XLim)*pos(3) + pos(1);
y2norm = @(y)(y-ax.YLim(1))/diff(ax.YLim)*pos(4) + pos(2);

hold on; 
annotation('doublearrow',x2norm([tHM1 tHM2]),y2norm([t/2 t/2]), ...
    'head1length',5,'head2length',5,'linewidth',2)
text(tHM1-0.05,t/2,'FWHM','fontweight','bold','fontsize',13, ... 
    'horizontalalignment','right')

annotation('doublearrow',x2norm([t3M1 t3M2]),y2norm([2*t/3 2*t/3]), ...
    'head1length',5,'head2length',5,'linewidth',2)
text(t3M1-0.05,2*t/3,'FW3M','fontweight','bold','fontsize',13, ...
    'horizontalalignment','right')

plot(tG1,tan1,'k','linewidth',3); scatter((tTro+2)*k,egAP(tTro+2),'r','filled')
text(tG1(end),egAP(tTro+2),'dV(0.06)','fontweight','bold','color','r','fontsize',13)
plot(tG2,tan2,'k','linewidth',3); scatter((tTro+12)*k,egAP(tTro+12),'r','filled')
text(tG2(2),egAP(tTro+12),'dV(0.4)','fontweight','bold','color','r','fontsize',13, ... 
    'horizontalalignment','right')

plot([tTro tTro]*k,[ax.YLim(1) t],'--','color','k','linewidth',1)
plot([tPeak tPeak]*k,[ax.YLim(1) p],'--','color','k','linewidth',1)
annotation('doublearrow',x2norm([tTro tPeak]*k),y2norm([min(ylim) min(ylim)]))
text([tTro+tPeak]*k/2,ax.YLim(1),'TtP','fontweight','bold','fontsize',13, ... 
    'HorizontalAlignment','center','verticalalignment','bottom')

title(sprintf('Neuron %d from dset %d',  ... 
    feets.whichCell(whichAP,1),feets.whichCell(whichAP,2)))
hold off; axis off;