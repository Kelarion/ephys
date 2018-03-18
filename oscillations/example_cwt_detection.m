%% load the data and LFP
if ~exist('lief','var')
    dataFolder = 'C:\DATA\Spikes\Janelia\muad_dib\2140_wheel_6-19-17\';
    sp = loadJRCdir(dataFolder);
    dat = load_channel(sp.dtype,[dataFolder sp.dat_path],sp.n_channels_dat, ... 
        sp.chanMap(sp.mainChannel(3)),1,sp.nSampDat);
    lief = getLFP(dat,sp.sample_rate);
end

[WT, F, newt, COI] = cwtnarrow(lief,sp.sample_rate,[50 15],'voicesperoctave',20);
F = F(:);
 
tdat = [1:length(dat)]/sp.sample_rate;

%% plot
magthresh = 150; % threshold on the magnitude
freqthresh = 2; % minimum number of frequency bands
isev = sum(abs(WT)>=magthresh,1)>=freqthresh;
timbs = findRipples(WT,F,T);

figure('units','normalized','position',[0.0406 0.3269 0.9406 0.5778]);
ax(1) = subtightplot(2,1,1,[0.005,0.01]);
ax(1).YLim = [min(dat) max(dat)];
shadeEvents(ax(1),timbs')
hold on
plot(ax(1),tdat(1:2:end),dat(1:2:end))
hold off
ax(1).XTick = [];
ax(1).YTick = [];

ax(2) = subtightplot(2,1,2,[0.005,0.01]);
imagesc(newt,[],abs(WT))
yt = yticks; % fix labels
flabs = split(num2str(F(yt)',3));
ax(2).YTickLabels = flabs;
xlabel('Time (sec)')

% ax(3) = subtightplot(3,1,3,[0.005,0.01]);
% stem(newt(evOn),width)

ax(1).XLim = [min(tdat) max(tdat)];
ax(2).XLim = [min(tdat) max(tdat)];
% ax(3).XLim = [min(tdat) max(tdat)];
linkaxes(ax,'x');

%%

width = zeros(size(evOn));
for ii = 1:nEvs
    thisEv = mag(:,evOn(ii):evOff(ii));
    [~,indf] = max(max(thisEv,[],2)); % which frequency was the peak
    [~,indt] = max(max(thisEv,[],1));
    m = max(thisEv(:,indt));
    [~,i1]  = min(abs(thisEv(1:indf,indt)-(m/2)));
    [~,i2]  = min(abs(thisEv(indf:end,indt)-(m/2)));
    i2 = i2+indf;
    width(ii) = i2 - i1;
end