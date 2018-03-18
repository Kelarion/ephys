%% load the data and LFP
if ~exist('lief','var')
    dataFolder = 'C:\DATA\Spikes\Janelia\muad_dib\2140_wheel_6-19-17\';
    sp = loadJRCdir(dataFolder);
    dat = load_channel(sp.dtype,[dataFolder sp.dat_path],sp.n_channels_dat, ... 
        sp.chanMap(sp.mainChannel(3)),1,sp.nSampDat);
    lief = getLFP(dat,sp.sample_rate);
end
slice = 1144900:1279900; % change this as you want

%% set wavelet parameters
fb = [2; 4; 8]; % variance of gaussian envelope, in seconds
k0 = 25; % center frequency of sinusoid (Hz)

%% iterate through parameters and plot
wt = zeros(length(slice),length(fb),length(k0));
win = zeros(length(fb),2);

fig = figure('Units','Normalized','Position',[0.0953 0.0417 0.6781 0.8796]);
tim = [1:length(slice)]/sp.sample_rate;
for ja = 1:length(fb)
    win(ja,:) = [-1 1]*(k0);
    nsamp = 2*sp.sample_rate;
    [psi, x] = cmorwavf(win(ja,1),win(ja,2),nsamp,fb(ja),1); % gabor wavelet
    x = x/k0;
    psi = psi/sum(abs(psi));
    wt(:,ja) = conv(dat(slice),conj(psi),'same');
    
    % plot
    powax(ja) = subplot(length(fb),2,2*ja-1);
    plot(tim,abs(wt(:,ja)),'linewidth',2)
    hold on
    plot(tim,real(wt(:,ja)),'linewidth',0.5)
    winsize = diff(win(ja,:));
    title(['Power using ' num2str(fb(ja)) ' variance'])
    if ja<length(fb), powax(ja).XTick = []; end
    
    wavax(ja) = subplot(length(fb),2,2*ja);
    plot(x,abs(psi),'linewidth',2)
    hold on; plot(x,real(psi),'color','r')
    plot(x,imag(psi),'--','color','r')
    title('Corresponding wavelet')
    wavax(ja).XLim = [-max(fb)/k0 max(fb)/k0];
    if ja<length(fb), wavax(ja).XTick = []; end
end
linkaxes(powax,'x')

% compare to what the MatLab function gives you
% [matwt, f] = cwt(lief(slice), 'amor');









