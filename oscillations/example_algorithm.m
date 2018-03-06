% example
dur = 120; %seconds

JFolder = 'C:\DATA\Spikes\muad_dib\2217_wheel_6-19-17\';
sp = loadJRCdir(JFolder,true);
goodcells = contains(sp.notes,'single');
goodchans = sp.mainChannel(goodcells);
egchan = goodchans(2);

c = firpmord([15 20 40 45], [0 1 0], [0.01 0.001 0.01], ...
    sp.sample_rate, 'cell');
if mod(c{1},2) == 1, c{1} = c{1} + 1; end
filt_low = firpm(c{:});

%% compute and plot
raw_dat = load_channel(sp.dtype,[JFolder sp.dat_path], ... 
    sp.n_channels_dat,sp.chanMap(egchan),1,sp.sample_rate*dur);
laschan = load_channel(sp.dtype,[JFolder sp.dat_path], ... 
    sp.n_channels_dat,65,1,sp.sample_rate*dur);
opto = cleanAux(laschan);

lief = getLFP(raw_dat,sp.sample_rate,filt_low);
% d = designfilt('highpassfir', ...    % Response type
%     'StopbandFrequency',500, ...     % Frequency constraints
%     'PassbandFrequency',650, ...
%     'StopbandAttenuation',55, ...    % Magnitude constraints
%     'PassbandRipple',4, ...
%     'DesignMethod','kaiserwin', ...  % Design method
%     'ScalePassband',false, ...       % Design method options
%     'SampleRate',sp.sample_rate)  ;
% ape = fftfilt(d,raw_dat);

f1 = 31; %Hz
K1 = sin(2*f1*linspace(0,pi*0.25,25000*0.25));
K1 = K1.*gausswin(25000*0.25)'/sum(gausswin(25000*0.25));
dinger1 = conv(lief,K1,'same');
hinger1 = hilbert(dinger1);

f2 = 25; %Hz
K2 = sin(2*f2*linspace(0,pi*0.25,25000*0.25));
K2 = K2.*gausswin(25000*0.25)'/sum(gausswin(25000*0.25));
dinger2 = conv(lief,K2,'same');
hinger2 = hilbert(dinger2);

%% plot
fig = figure('Units', 'Normalized','Position', [0.0484 0.3028 0.8042 0.6019]);

ax1 = subplot(3,1,1); 
ax1.Position = [0.1300 0.6722 0.7750 0.2800];
ax1.YLim = [min(raw_dat) max(raw_dat)];
shadeEvents(ax1,1:length(opto),opto,[0.8 0.8 0.8], 'FaceAlpha',0.8,'linestyle','none')
hold on; plot(raw_dat)
% plot(lief*2,'linewidth',2)
axis off

ax2 = subplot(3,1,2);
ax2.Position = [0.1300 0.3572 0.7750 0.3042];
plot(dinger2)
hold on; plot(abs(hinger2),'linewidth',2)
axis off

ax3 = subplot(3,1,3);
ax3.Position = [0.1300 0.0628 0.7750 0.2818];
plot(dinger1)
hold on; plot(abs(hinger1),'linewidth',2)
axis off

linkaxes([ax1 ax2 ax3],'x')
ax1.XLim = [522990   650920];
% ax1.YLim = [-5126 4168];

