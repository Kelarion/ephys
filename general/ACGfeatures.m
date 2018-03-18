function af = ACGfeatures(ACG,ACG_bins)
% af = ACGfeatures(ACG,ACG_bins)
%
% Inputs:
%  - ACG [nCell,nBin]: Autocorrelograms of all cells
%  - ACG_bins [nCell or 1,nBin]: Time bins of corresponding ACGs
% Output struct fields:
%   > t_ref: Refractory period, measured as time to half-maximum
%   > peakiness: Maximum/asymptotic value
%   > PW2: Width of peak at half-maximum
%   > PW3: width of peak at 2/3-maximum
% To add:
%   > [p_max, f_max]: Power and frequency of strongest Fourier component

nClu = size(ACG,1);
upsample_peak = 2;
upsample_risetime = 5;

rfp = zeros(nClu,1); % refractoriness
r_max = zeros(nClu,1); % maximum firing rate
r_inf = zeros(nClu,1); % asymptotic firing rate
r_mean = zeros(nClu,1);
t_mean = zeros(nClu,1); % mean lag
PW2 = nan(nClu,1); % peak width at half-maximum
PW3 = nan(nClu,1); % peak width at 2/3-maximum
for iClu = 1:nClu 
    acg = ACG(iClu,:); % normalise
    if size(ACG_bins,1) == 1 
        t = ACG_bins;
    else
        t = ACG_bins(iClu,:);
    end
    Fs = 1/mean(diff(t));
    th = round(length(t)/2);
    t100 = t <= 0.1;
    
    acg_avg = smooth(acg,10);
    acg_smooth = smooth(acg,0.02,'loess');
    
    t_mean(iClu) = dot(acg(t100),t(t100))/sum(acg(t100));
    
    [~,tpk] = max(acg(1:th));
    if tpk == 1 
%         warning(['Cell ' num2str(iClu) ' not refractory, skipping'])
        rfp(iClu) = NaN; 
        r_max(iClu) = NaN; 
        r_inf(iClu) = NaN;
        r_mean(iClu) = NaN;
        t_mean(iClu) = NaN;
        continue
    end
    if tpk >= th/2
        r_max(iClu) = mean(acg_avg(1:th)); % if the peak is too late, assume it's not real
    else
        r_max(iClu) = mean(acg(tpk-1:tpk+1));
    end
    r_inf(iClu) = mean(acg_avg(end-50:end));
    r_mean(iClu) = mean(acg_avg);
    
    if r_max(iClu) > 1.8*r_inf(iClu)
        tq = interp1(1:length(t),t,1:(1/upsample_peak):th);
        d = find(tq == t(tpk));
        acg_interp = interp1(t,acg_avg,tq);
        [~,t1] = min(abs(acg_interp(1:d) - r_max(iClu)/2));
        [~,t2] = min(abs(acg_interp(d:end) - r_max(iClu)/2));
        PW2(iClu) = (t2+d - t1)/(upsample_peak*Fs);
        
        [~,t1] = min(abs(acg_interp(1:d) - 2*r_max(iClu)/3));
        [~,t2] = min(abs(acg_interp(d:end) - 2*r_max(iClu)/3));
        PW3(iClu) = (t2+d - t1)/(upsample_peak*Fs);
    end
    
    tq = interp1(1:length(t),t,1:(1/upsample_risetime):tpk);
    acg_interp = interp1(t,acg,tq);
    thm = find(acg_interp >= 0.85*r_max(iClu),1,'first');
    rfp(iClu) = thm/(upsample_risetime*Fs);
    
    
end

af.t_ref = rfp;     af.t_mean = t_mean;
af.peakiness = r_max./r_inf;
af.r_mean = r_mean; 
af.PW2 = PW2;       af.PW3 = PW3;

end

