function sf = spikeFeatures(wf,Fs)
% function sf = spikeFeatures(wf,Fs)
% 
% Compute some features of the extracellular action potential (EAP)
% waveform, which are likely to relate to the density of different ion
% conductances (Gold, alia, Buszaki, 2006). Note that spikes for which
% these features can't easily be computed (i.e. badly-shaped spikes) have
% NaN values instead.
% 
% Inputs
% - wf [nSamp,nClu(,nSpikes)]: The EAP waveforms.
% (- wf [nSamp,nTrode,nSpikes]: alternative if computing on single spikes.)
% - Fs [1,1]: Sample rate (default 1).
%
% Output struct
% - sf.TP [nClu,1]: The mean trough-to-peak time and variance.
% - sf.FWXM [nClu,2]: The full-width at x-maximum, a measure of the width
%   of the sodium current. Column 1 is half-max, column two is 1/3 max.
% - sf.Grad [nClu,2]: The gradient of the waveform at 0.067 and 0.4 ms
%   after the trough time. The specific times can be changed in function.
% - sf.Ratio [nClu,1]: The ratio of minimum to maximum deflection from
%   baseline, i.e. abs(peak)/abs(trough).
% - sf.Amps [nClu,1]: The peak-trough values.
% - sf.I_cap [nClu,1]: Amplitude of the capacitative phase relative to the 
%   trough, which is NaN if there is no prominent capacitative phase. But
%   this most likely depends on electrode position.
%
% Features are based on discussion with I-Chun Lin, Cortex Lab UCL

useRaw = size(wf,3) > 1;
nClu = size(wf,2);
if useRaw 
    nSpk = size(wf,3);
else 
    nSpk = 1;
end
nSamp = size(wf,1);

if nargin<2, Fs = 1; end

tGrad = [0.067 0.6]; % decide where to take the gradient

Grad   = nan(nClu, 2, nSpk);   % gradient @0.067 and 0.4 ms
TP    = nan(nClu, 1, nSpk);   % trough-to-peak time
Ratio  = nan(nClu, 1, nSpk); % Height ratio (peak/trough)
Amps  = nan(nClu, 1, nSpk); % Height ratio (peak/trough)
ICap = nan(nClu, 1, nSpk);    % Capatitative current size (cap/trough)
FWHM   = nan(nClu, 1, nSpk); % width at 1/2 maximum
FW3M   = nan(nClu, 1, nSpk); % width at 1/3 maximum

for iClu = 1:nClu
    if useRaw 
        % important: I assume that the spikes are smoothed beforehand!
        theseAP = permute(wf(:,iClu,:),[1,3,2]); 
        v0 = theseAP(1,:); % baseline voltage, is 0 for templates
        theseAP = theseAP - v0;
        clear v0
        
        [p, tPeak] = max(theseAP,[],1);
        [t, tTro] = min(theseAP,[],1);
        
        % quality control
        toobig = abs(t) > quantile(abs(t),0.99);
        toosmall = abs(t) < quantile(abs(t),0.02);
        tooweird = (tTro < 5) | (tTro > 15);
        isSpike = ~(toobig | toosmall | tooweird);
        clear toobig toosmall
        
        isCapac = tPeak < tTro;  % if it has a strong capacitative phase
%         for ispk = find(isCapac & isSpike)
%             [p(ispk),tptemp] = max(theseAP(tTro(ispk):end) - v0(ispk));
%             tPeak(ispk) = tptemp + tTro(ispk) - 1;
%         end
        
        finalV0 = mean(theseAP(end-4:end,:));
        isBiphas = abs(p-finalV0) >= 0.2;
        goodSpike = (~isCapac) & isSpike;
%         goodSpike = isBiphas & isSpike;
        clear finalV0 isCapac isBiphas isSpike
        
        % FWHM/ FW3M
        % to efficiently compute FWHM, we assume a linear interpolation of
        % the waveform to find the exact half-max crossings
        halfmax = -abs(t(goodSpike))/2;
        [tH1,tH2,toowideH] = interpfind(halfmax,theseAP(:,goodSpike),Fs);
        
        thirdmax = -abs(2*t(goodSpike))/3;
        [tT1,tT2,toowideT] = interpfind(thirdmax,theseAP(:,goodSpike),Fs);
        
        inBoth = (~toowideH) & (~toowideT);
        keepSpikes = goodSpike; % combine logical vectors of different size
        keepSpikes(goodSpike) = inBoth;
        
        FWHM(iClu,1,keepSpikes) = tH2(~toowideT(~toowideH)) - tH1(~toowideT(~toowideH));
        FW3M(iClu,1,keepSpikes) = tT2(~toowideH(~toowideT)) - tT1(~toowideH(~toowideT));
        
        % Gradients
        ngood = sum(keepSpikes);
        d1 = round(tGrad(1)*Fs/1000); d2 = round(tGrad(2)*Fs/1000);
        dv = [zeros(1,ngood); diff(theseAP(:,keepSpikes))];
        ind1 = sparse(tTro(keepSpikes)+d1,1:ngood,true(1,ngood),nSamp,ngood);
        ind2 = sparse(tTro(keepSpikes)+d2,1:ngood,true(1,ngood),nSamp,ngood);
        Grad(iClu,1,keepSpikes) = dv(ind1)/(Fs/1000); 
        Grad(iClu,2,keepSpikes) = dv(ind2)/(Fs/1000);
        
        % Other things
        TP(iClu,1,keepSpikes) = (tPeak(keepSpikes) - tTro(keepSpikes))/Fs;
        Ratio(iClu,1,keepSpikes) = abs(p(keepSpikes))./abs(t(keepSpikes));
        Amps(iClu,1,keepSpikes) = p(keepSpikes) - t(keepSpikes);
        
    else
        thisAP = smooth(wf(:,iClu),0.1,'loess');
        v0 = median(wf(1:10,iClu)); % baseline voltage, is 0 for templates
        [p, tPeak] = max(thisAP - v0);
        [t, tTro] = min(thisAP - v0);
        if tPeak < tTro  % if it has a strong capacitative phase
            c = p;
            ICap(iClu) = abs(c)/abs(t);
            [p,tptemp] = max(thisAP(tTro:end) - v0); % we want the second peak
            tPeak = tptemp + tTro - 1;
        end 
        finalV0 = mean(thisAP(end-4:end));
        if abs(p-finalV0) >= 0.2 % threshold whether positive phase is large enough
            TP(iClu) = (tPeak - tTro)/Fs; % trough-to-peak (sec)
            Ratio(iClu) = abs(p)/abs(t);   % height-ratio
        end
        
        t_temp = tTro-5:tTro+9; tq = tTro -5:0.1:tTro+9;
        AP_interp = interp1(t_temp,thisAP(t_temp),tq);
        d = find(tq == tTro);
        [~, H1] = min(abs(AP_interp(1:d)+abs(t)/2));
        [~, H2] = min(abs(AP_interp(d:end)+abs(t)/2));
        tH1 = tq(H1);tH2 = tq(H2 +d -1); % FWHM
        FWHM(iClu) = (tH2-tH1)/Fs;
        
        [~, H1] = min(abs(AP_interp(1:d)+abs(2*t)/3));
        [~, H2] = min(abs(AP_interp(d:end)+abs(2*t)/3));
        tH1 = tq(H1);tH2 = tq(H2 +d -1); % FW3M
        FW3M(iClu) = (tH2-tH1)/Fs;
        
        t1 = round(tGrad(1)*Fs/1000); t2 = round(tGrad(2)*Fs/1000); % grad points
        dv = [0; diff(thisAP)];
        Grad(iClu,1) = dv(tTro+t1); Grad(iClu,2) = dv(tTro+t2);
        
    end
    
end

sf.TP = TP;
sf.FWXM = [FWHM FW3M];
sf.Grad = Grad;
sf.Ratio = Ratio;
sf.Amps = Amps;
sf.I_cap = ICap;

%% helper function to find times of half-max and 1/3-max
function [t1,t2,biguns] = interpfind(vtarget,spikes,Fs)
    % [t1,t2,grads,biguns] = interpfind(vtarget,spikes,Fs)
    %
    % Find t s.t. spikes(t(i),i) = vtarget(i), using a linear
    % interpolation of `spikes`. Assumes only two crossings of vtarget,
    % used to find times of half-max and 1/3-max crossings.
    
    [nt,~] = size(spikes);
    
    [ro,co] = find(spikes <= vtarget);
    [~,ic] = unique(co); % get the first non-zero entry for each spike
    a1 = ro(ic) - 1; a2 = a1 + 1; % points on either side of first crossing
    b1 = ro([ic(2:end)-1;end]); b2 = b1 + 1; % and of second crossing
    
    biguns = (b1 - a1) > Fs*(1e-3); % no things wider than 1 ms
    a1 = a1(~biguns); a2 = a2(~biguns);
    b1 = b1(~biguns); b2 = b2(~biguns);
    vtarget = vtarget(~biguns);
    S = spikes(:,~biguns);
    n = sum(~biguns); 
    
    x = 1:n;
    va = reshape(S(sparse([a1;a2]', [x,x], true(2*n,1),nt,n)),2,n);
    vb = reshape(S(sparse([b1;b2]', [x,x], true(2*n,1),nt,n)),2,n);
    dtdv_a = -1./diff(va); % t1 = t_a1 + (v_a1 - v_hm)(dt/dv)_a1
    dtdv_b = -1./diff(vb); % similar as above
    t1 = a1' + (va(1,:) - vtarget).*dtdv_a;
    t2 = b1' + (vb(1,:) - vtarget).*dtdv_b;
end

end