%% set filter
c = firpmord([15 20 40 45], [0 1 0], [0.01 0.001 0.01], 25000, 'cell');
if mod(c{1},2) == 1, c{1} = c{1} + 1; end
filt_low = firpm(c{:});

%% load the data and LFP
dataFolder = 'C:\DATA\Spikes\Janelia\muad_dib\2217_wheel_6-19-17\';
sp = loadJRCdir(dataFolder);
dat = load_channel(sp.dtype,[dataFolder sp.dat_path],sp.n_channels_dat, ...
    sp.chanMap(sp.mainChannel(3)),1,sp.nSampDat);
lief = getLFP(dat,sp.sample_rate,filt_low);

[WT, F, newt, COI] = cwtnarrow(lief,sp.sample_rate,[50 15],'voicesperoctave',20);
F = F(:);
 
tdat = [1:length(dat)]/sp.sample_rate;

%% plotdat 
timbs = findRipples(WT,F,newt);

figure('units','normalized','position',[0.0406 0.3269 0.9406 0.5778]);
ax(1) = subtightplot(2,1,1,[0.005,0.01]);
ax(1).YLim = [min(dat) max(dat)];
shadeEvents(timbs',[],ax(1))
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
