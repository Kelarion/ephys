%% set filter
c = firpmord([15 20 40 45], [0 1 0], [0.01 0.001 0.01], 25000, 'cell');
if mod(c{1},2) == 1, c{1} = c{1} + 1; end
filt_low = firpm(c{:});

%% load the data and LFP
dataFolder = 'C:\DATA\Spikes\Janelia\muad_dib\2140_wheel_6-19-17\';
sp = loadJRCdir(dataFolder);
dat = load_channel(sp.dtype,[dataFolder sp.dat_path],sp.n_channels_dat, ...
    sp.chanMap(sp.mainChannel(3)),1,sp.nSampDat);
lief = getLFP(dat,sp.sample_rate,filt_low);

[WT, F, newt, COI] = cwtnarrow(lief,sp.sample_rate,[50 15],'voicesperoctave',30);
F = F(:);

berty = hilbert(lief);
phas = angle(berty);
mag = abs(berty);

[timbs,freqs,pkmag] = findRipples(WT,F,newt,params);
%% plot
tdat = [1:size(rawDat.Data.x,2)]/sp.sample_rate;
tlfp = [1:length(lief)]/sp.sample_rate;

[~,ii] = sort(freqs);
inds = findNearestPoint(sort(freqs),[sort(F); nan(1000,1)]); % because this function is weird
freqinds = zeros(size(freqs));
freqinds(ii) = length(F) - inds;
freqtims = mean(timbs,2);

figure('units','normalized','position',[0.0406 0.3269 0.9406 0.5778]);
ax(1) = subtightplot(2,1,1,[0.005,0.01]);
ax(1).YLim = [min(rawDat.Data.x(chan,1:2:end)) max(rawDat.Data.x(chan,1:2:end))];
shadeEvents(timbs',[],ax(1))
hold on
plot(ax(1),tdat(1:2:end),rawDat.Data.x(chan,1:2:end))
% plot(ax(1),tlfp(1:2:end),phas(1:2:end).*mag(1:2:end))
hold off
ax(1).XTick = [];
ax(1).YTick = [];

ax(2) = subtightplot(2,1,2,[0.05,0.01],[0.1, 0.05],[0.05 0.05]);
imagesc(newt,[],abs(WT))
caxis([0 max(pkmag)])
yt = yticks; % fix labels
flabs = split(num2str(F(yt)',3));
ax(2).YTickLabels = flabs;
xlabel('Time (sec)')
ylabel('Frequency')
hold on;
text(freqtims,freqinds,split(num2str(freqs',3)), ...
    'HorizontalAlignment','center','fontweight','bold')

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

