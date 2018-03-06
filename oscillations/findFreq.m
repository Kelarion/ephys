function [rippleon, rippleoff] = findFreq(raw_channel, leaf, tfrq, thr, fb)
% [rippleon, rippleoff] = findFreq(raw_channel, leaf, tfrq, thr[,fb])
%
% raw_channel: raw AP-band signal
% leaf: LFP, narrow-band filtered near target frequency
% tfrq: target frequency
% thr: threshold on frequency magnitude (units of dB)
% minDur: minimum duration of event (units of cycles) 
% fb: variance of gabor wavelet envelope



win = [-1 1]*tfrq;
nsamp = 2*sp.sample_rate;
[psi, x] = cmorwavf(win,win,nsamp,fb,1); % gabor wavelet
x = x/tfrq;
psi = psi/sum(abs(psi));

wt = conv(leaf,conj(psi),'same');
hinger = hilbert(wt);

putative = abs(hinger) >= thr*(10^5);
rippleon = find(diff(putative) > 0);
rippleoff = find(diff(putative) < 0);
if rippleon(1) > rippleoff(1), rippleoff(1) = []; end
if length(rippleon) > length(rippleoff), rippleon(end) = []; end

if nargout < 1 
    fig = figure('Units', 'Normalized','Position', [0.0484 0.3028 0.8042 0.6019]);
    ax1 = subplot(2,1,1);
    ax2 = subplot(2,1,2);
    ax1.XTick = [];
    whitebg(fig,'k');
    goodtotal = 0;
    h = hilbert(leaf);
    t = 1:length(abs(hinger));
    dv = abs(diff(raw_channel,2));
    for i = 1:length(rippleon)
        istart = rippleon(i);
        ifin = rippleoff(i);
        tea = t(istart:ifin);
        rah = raw_channel(istart:ifin);
        %     truspks = logical(conv(logical(all_spks{3}(4,istart:ifin)),ones(1, 40),'same'));
        
        if ifin - istart >= 3906 && ~any(dv(istart:ifin) > 2500)
            % must be more than 4 periods
            ah = max(leaf(istart:ifin))/(2*pi);
            
            plot(ax1,tea,rah)
            plot(ax2,t(istart:ifin),leaf(istart:ifin))
            %         hold(ax1,'on'); plot(ax1,tea(truspks),rah(truspks)); hold(ax1,'off');
            hold(ax2,'on'); plot(ax2,t(istart:ifin), angle(h(istart:ifin))*ah); hold(ax2,'off');
            linkaxes([ax1 ax2], 'x')
            goodtotal = goodtotal + 1;
            pause
        end
    end
end

end
