function sf = spikeFeatures(wf,Fs)
% function sf = spikeFeatures(wf,Fs)
% 
% Compute some features of the extracellular action potential (EAP)
% waveform, which are likely to relate to the density of different ion
% conductances (Gold, alia, Buszaki, 2006). If this is computed on the
% raw waveforms, then the variance of those features is also returned; if
% they are computed on the mean waveform or template (from spike-sorting) 
% then just the quantities are returned. (At the moment it only works for
% the latter scenario.)
% 
% Inputs
% - wf [nSamp,nClu(,nSpikes)]: The EAP waveforms.
% - Fs [1,1]: Sample rate (default 1).
%
% Output struct
% - sf.TP [nClu,1]: The mean trough-to-peak time and variance.
% - sf.FWXM [nClu,2]: The full-width at x-maximum, a measure of the width
%   of the sodium current. Column 1 is half-max, column two is 1/3 max.
% - sf.Grad [nClu,2]: The gradient of the waveform at 0.067 and 0.4 ms
%   after the trough time. The specific times can be changed in function.
% - sf.Amps [nClu,1]: The ratio of minimum to maximum deflection from
%   baseline, i.e. abs(peak)/abs(trough).
% - sf.I_cap [nClu,1]: Amplitude of the capacitative phase relative to the 
%   trough, which is NaN if there is no prominent capacitative phase. But
%   this most likely depends on electrode position.
% - sf.'field'_var (optional): The variance associated with all of the
%   measures above, can take substantially longer to compute.
%
% Feature computing code contributed by I-Chun Lin, Cortex Lab UCL


useRaw = size(wf,3) > 1;
nClu = size(wf,2);
if useRaw 
    nSpk = size(wf,3);
end

if nargin<2, Fs = 1; end

tGrad = [0.067 0.4]; % decide where to take the gradient

Grad   = nan(nClu, 2);   % gradient @0.067 and 0.4 ms
TP    = nan(nClu, 1);   % trough-to-peak time
Ratio  = zeros(nClu, 1); % Height ratio (peak/trough)
ICap = nan(nClu, 1);    % Capatitative current size (cap/trough)
FWHM   = zeros(nClu, 1); % width at 1/2 maximum
FW3M   = zeros(nClu, 1); % width at 1/3 maximum

for iClu = 1:nClu
    if useRaw % don't know if I'll ever use this
        
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
sf.Amps = Ratio;
sf.I_cap = ICap;

