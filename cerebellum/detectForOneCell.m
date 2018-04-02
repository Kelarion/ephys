ksRoot = 'D:\data\';

params.pThresh = 170; % arbitrary magnitude units
params.tThresh = 4; % periods
params.minBetween = 0.5; % seconds
params.noiseCutoff = 18; % Hz
nvox = 20; % number of voices

if ~exist('filt_low','var') % takes a long time to make this filter
    c = firpmord([15 20 40 45], [0 1 0], [0.01 0.001 0.01], 25000, 'cell');
    if mod(c{1},2) == 1, c{1} = c{1} + 1; end
    filt_low = firpm(c{:}); 
end % but it seems to give the best performance for such a narrow bandwidth

cid = 4;
d = 3;
k = 3;
%% load
parentDir = [ksRoot '\' db(k).name '\'];
subDir = dir([parentDir '*' db(k).date]);

if length(subDir) > 1
    subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
    if isempty(subDir)
        disp('can''t find it ¯\_(?)_/¯')
    end
end
dataDir = [parentDir subDir.name '\'];

sp = loadJRCdir(dataDir);
rawDat = memmapfile([dataDir sp.dat_path], 'Format',  ...
    {sp.dtype, [sp.n_channels_dat sp.nSampDat], 'x'});

%% detect
if db(k).hasOpto(d)
    if db(k).badAux
        opto = cleanAux(rawDat.Data.x(66,:)); % get laser times
    else
        op = rawDat.Data.x(66,:)>0.5*max(rawDat.Data.x(66,:));
        t = 1:length(op);
        optoon = find(diff([0 op])>0);
        optooff = find(diff([0 op])<0);
        
        opto = [optoon(:) optooff(:)];
        
        tooSoon = find(diff(opto,[],2) < 50); % events occuring too soon after the last
        opto(tooSoon-1,2) = opto(tooSoon,2); % merge with adjacent events
        opto(tooSoon,:) = [];
        opto = t(opto);
    end
    optoTimes = opto/sp.sample_rate;
    dopto = diff(optoTimes(:,1));
    inds = find(dopto > 1.2*min(dopto));
    boutStart = [optoTimes(1,1); optoTimes(inds+1,1)];
    boutEnd = [optoTimes(inds,2); optoTimes(end,2)];
    stim = [boutStart boutEnd]; % the 'bouts' of stimulation
    
    optoFreq = zeros(length(stim),1);
    wr = WithinRanges(optoTimes(:,1),stim,1:length(stim),'vector');
    optoFreq(wr,:) = 1./[dopto(1); dopto]; % freq of each bout
end

chan = sp.chanMap(sp.mainChannel(sp.cids == cid));
lief = getLFP(rawDat.Data.x(chan,:),sp.sample_rate,filt_low);
phas = angle(hilbert(lief));

[WT, F, newt, COI] = cwtnarrow(lief,sp.sample_rate,[50 15],'voicesperoctave',nvox);
F = F(:);
[timbs,freqs,pkmag] = findRipples(WT,F,newt,params);

%% compute
if isempty(timbs), error('no events'); end

if db(k).hasOpto(d)
    if size(timbs,1) > size(stim,1)
        inds = findNearestPoint(mean(stim,2),mean(timbs,2));
        detStims = 1:size(stim,1);
    else
        couldveFound = optoFreq > 2;
        detStims = unique(findNearestPoint(mean(timbs,2),mean(stim,2)));
        inds = findNearestPoint(mean(stim(detStims,:),2),mean(timbs,2));
    end
    
    % find which were induced by laser:
    % consider events which overlap with stimulation as induced
    overlap = timbs(inds,1)<stim(detStims,2) | timbs(inds,2)>stim(detStims,1);
    isStim = ismember(1:size(timbs,1),inds(overlap));
    isStim = isStim(:);
    sfreq = nan(length(isStim),1); % frequency of stim when it caused a ripple
    sfreq(inds,:) = optoFreq(detStims(overlap));

else
    isStim = false(size(timbs,1),1);
    sfreq = nan(size(timbs,1),1);
end
evInds = [1:size(timbs,1)]';

stInds = round(sp.st(sp.clu == cid)*sp.sample_rate);
wr = WithinRanges(sp.st(sp.clu == cid),timbs,1:size(timbs,1),'vector');
p = phas(stInds(logical(wr))); % get the phase of each spike relative to the events

