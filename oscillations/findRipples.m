function [evTimes,peakFreq,peakPower] = findRipples(WT,F,T,params)
% [evTime,peakFreq,peakPower] = findRipples(WT,F,T[,params])
%
% WT (nF,nT) is wavelet transform, or transform magnitude 
% F and T are the frequency and time labels
% 
% params has 4 accepted fields:
% - pThresh is threshold on magnitude, i.e. |WT|^2 (default 150)
% - tThresh is minimum event duration in periods (default 3)
% - minBetween is minimum distance for two events to be considered different,
%   in the same units as T (default 0.5 seconds)
% - noiseCutoff is the minimum spread in frequency domain for an event to
%   be considered noise, in Hz (default 18)
%
% see 'cwt' or 'cwtnarrow', which generate the first three inputs

if nargin < 4, params = struct; end
if ~isfield(params,'pThresh'), params.pThresh = 150; end
if ~isfield(params,'tThresh'), params.tThresh = 3; end
if ~isfield(params,'minBetween'), params.minBetween = 0.5; end % assume seconds
if ~isfield(params,'noiseCutoff'), params.noiseCutoff = 18; end

magThr = params.pThresh;
minDur = params.tThresh;
minBetween = params.minBetween;
noiseThr = params.noiseCutoff;

dt = mean(diff(T));
Fs = 1/dt;

minSampBetween = ceil(minBetween*Fs);

mag = abs(WT);%/max(max(abs(WT)));

isev = sum(mag>=magThr,1)>=2;
evOn = find(diff([0 isev]) > 0);
evOff = find(diff([isev 0]) < 0);
nEvs = length(evOn);

peakPower = nan(nEvs,1);
peakFreq = nan(nEvs,1);
for iEv = 1:nEvs % go through each event and make sure it passes all criteria
    if isnan(evOn(iEv))
        continue
    end
    thisEv = mag(:,evOn(iEv):evOff(iEv));
    m = max(max(thisEv,[],2));
    [~,indt] = max(max(thisEv,[],1)); 

    if (fspread(thisEv(:,indt),F) > noiseThr) % artifacts are wider in frequency domain
        evOn(iEv) = nan; evOff(iEv) = nan;
        continue
    end
    
    [~,indf] = max(mean(thisEv,2));
    fmax = F(indf); % if two events of the same frequency are too close together
    if iEv < nEvs  % then we consider them the same event
        [~,i2] = max(mean(mag(:,evOn(iEv+1):evOff(iEv+1)),2));
        if (abs(F(i2)-fmax) <= 3) && (diff(evOn(iEv:iEv+1)) <= minSampBetween)
            evOff(iEv) = nan; evOn(iEv+1) = nan;
            peakFreq(iEv) = fmax;
            peakPower(iEv) = m;
            continue % the next event will be skipped
        end
    end
    
    minInds = minDur/(fmax*dt); % finally, remove if too brief
    if (evOff(iEv) - evOn(iEv) < minInds)
        evOn(iEv) = nan; evOff(iEv) = nan;
        continue
    end
    
    peakFreq(iEv) = fmax;
    peakPower(iEv) = m;
end

evTimes = [evOn(~isnan(evOn)); evOff(~isnan(evOff))]./Fs;
evTimes = evTimes';
peakFreq = peakFreq(~isnan(peakFreq));
peakPower = peakPower(~isnan(peakPower));

end


%% ---------------------------------------------------------------------
function width = fspread(x,f)
% width = fspread(x,f)
% Spread of x in terms of frequency

[m,ind] = max(x);

[~,i1]  = min(abs(x(1:ind)-(m/2)));
[~,i2]  = min(abs(x(ind:end)-(m/2)));
i2 = i2+ind;

width = f(i1) - f(i2);

end
