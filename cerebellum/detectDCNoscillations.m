% set parameters
ksRoot = 'C:\DATA\Spikes\Janelia\';

verbose = 1;

params.pThresh = 130; % arbitrary magnitude units
params.tThresh = 3; % periods
params.minBetween = 0.1; % seconds
params.noiseCutoff = 18; % Hz
nvox = 20; % number of voices

clear db;
wheel_db

if ~exist('filt_low','var') % takes a long time to make this filter
    c = firpmord([15 20 40 45], [0 1 0], [0.01 0.001 0.01], 25000, 'cell');
    if mod(c{1},2) == 1, c{1} = c{1} + 1; end
    filt_low = firpm(c{:}); 
end % but it seems to give the best performance for such a narrow bandwidth

%% run
R = struct;
whichCell = [];
whichEvent = [];
spkPhase = [];
spkMag = [];
isOpto = [];
evFreq = [];
evTime = [];
evMag = [];
stimFreq = [];

allEvTimes = [];
allEvFreqs = [];
allOptoFreq = [];
oscID = [];
nMiss = 0;
meanErr = nan(length([db.depth]),1);
j = 1;
tic
for k = 1:length(db)
    for d = 1:length(db(k).depth)
        %% load
        parentDir = [ksRoot '\' db(k).name '\'];
        subDir = dir([parentDir '*' db(k).date]);
        
        if verbose
            nicename = join(split(db(k).name,'_'));
            disp(['Now analyzing: ' nicename{:} ', ' db(k).depth{d} 'um ...']); 
        end
        if length(subDir) > 1
            subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
            if isempty(subDir)
                disp('can''t find it �\_(?)_/�')
%                 continue
            end
        end
        dataDir = [parentDir subDir.name '\'];
        
        sp = loadJRCdir(dataDir);
        rawDat = memmapfile([dataDir sp.dat_path], 'Format',  ...
            {sp.dtype, [sp.n_channels_dat sp.nSampDat], 'x'});
        
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
        
        isgood = cellfun(@ischar,sp.notes); % because sometimes I was stupid
        isgood(isgood) = contains(sp.notes(isgood),'single');
        nCells = sum(isgood);
        goodCids = sp.cids(isgood);

        tdat = [1:size(rawDat.Data.x,2)]/sp.sample_rate;
        
        %% detect
        if verbose, disp(['Finding events for ' num2str(nCells) ' channels ...']); end
        CIDs = [];
        myst = nan(1,nCells); % how many known events (opto-induced) were missed
        errs = nan(1,nCells); % mean time not detected by the algorithm
        for iClu = 1:nCells
            chan = sp.chanMap(sp.mainChannel(sp.cids == goodCids(iClu)));
            
            if isfield(sp,'lfp_path')
                lfpRaw = memmapfile([dataDir sp.lfp_path], 'Format',  ...
                    {sp.dtype, [sp.n_channels_dat sp.nSampLFP], 'x'});
                lief = getLFP(lfpRaw.Data.x(chan,:),sp.lfp_sample_rate,filt_low);
                Fs = sp.lfp_sample_rate;
            else
                lief = getLFP(rawDat.Data.x(chan,:),sp.sample_rate,filt_low);
                Fs = sp.sample_rate;
            end
            
            [WT, F, newt, COI] = cwtnarrow(lief,Fs,[50 15],'voicesperoctave',nvox);
            F = F(:);
            [timbs,freqs,pkmag] = findRipples(WT,F,newt,params);
            
            if isempty(timbs), continue; end
            
            berty = hilbert(lief);
            phas = angle(berty);
            mag = abs(berty);
            
            if db(k).hasOpto(d)
                if size(timbs,1) > size(stim,1)
                    inds = findNearestPoint(mean(stim,2),mean(timbs,2));
                    detStims = 1:size(stim,1);
                    myst(iClu) = 0;
                else
                    couldveFound = optoFreq > 2;
                    detStims = unique(findNearestPoint(mean(timbs,2),mean(stim,2)));
                    inds = findNearestPoint(mean(stim(detStims,:),2),mean(timbs,2));
                    myst(iClu) = (size(stim(couldveFound,:),1)-size(timbs,1));
                end
                
                % find which were induced by laser:
                % consider events which overlap with stimulation as induced
                overlap = (timbs(inds,1)<=stim(detStims,2)) & (timbs(inds,2)>=stim(detStims,1));
                isStim = ismember(1:size(timbs,1),inds(overlap));
                isStim = isStim(:); 
                sfreq = nan(length(isStim),1); % frequency of stim when it caused a ripple
                sfreq(inds(overlap),:) = optoFreq(detStims(overlap));
                
                myst(iClu) = myst(iClu) + sum(~overlap);
                err = abs(timbs(inds(overlap),1)-stim(detStims(overlap),1)) ...
                    + abs(timbs(inds(overlap),2)-stim(detStims(overlap),2));
%                 err = abs(timbs(inds(overlap),1)-stim(detStims(overlap),1));
                errs(iClu) = mean(err);
            else
                isStim = false(size(timbs,1),1);
                sfreq = nan(size(timbs,1),1);
                myst(iClu) = 0;
                errs(iClu) = 0;
            end
            evInds = [1:size(timbs,1)]';
            
            stInds = round(sp.st(sp.clu == goodCids(iClu))*Fs);
            wr = WithinRanges(sp.st(sp.clu == goodCids(iClu)),timbs,1:size(timbs,1),'vector');
            p = phas(stInds(logical(wr))); % get the phase of each spike relative to the events
            m = mag(stInds(logical(wr)));
            
            allEvTimes = [allEvTimes; timbs];
            allEvFreqs = [allEvFreqs; freqs];
            oscID = [oscID; goodCids(iClu)*ones(size(timbs,1),1) ...
                d*ones(size(timbs,1),1) k*ones(size(timbs,1),1)];
            
            spkPhase = [spkPhase; p(:)];
            spkMag = [spkMag; m(:)];
            isOpto = [isOpto; isStim(wr(wr>0),:)];
            evFreq = [evFreq; freqs(wr(wr>0),:)];
            evTime = [evTime; timbs(wr(wr>0),:)];
            evMag = [evMag; pkmag(wr(wr>0),:)];
            stimFreq = [stimFreq; sfreq(wr(wr>0),:)];
            allOptoFreq = [allOptoFreq; sfreq(:)];
            CIDs = [CIDs; goodCids(iClu)*ones(length(p),1)];
            whichEvent = [whichEvent; evInds(wr(wr>0),:)];
        end
        if verbose, disp(['     done in ' num2str(toc) ' seconds']); end
        whichCell = [whichCell; CIDs d*ones(length(CIDs),1) k*ones(length(CIDs),1)];
        nMiss = nMiss + nansum(myst);
        meanErr(j) = nanmean(errs);
        j = j+1;
    end
end

R.osc.Times = allEvTimes;
R.osc.Freq = allEvFreqs;
R.osc.optoFreq = allOptoFreq;
R.osc.whichCell = oscID;  

R.whichEvent = whichEvent;
R.spkPhase = spkPhase;  
R.spkMag = spkMag;
R.magnitude = evMag;
R.freq = evFreq;        
R.interval = evTime;
R.isOpto = logical(isOpto);      
R.optoFreq = stimFreq;  
R.datasets = db;        
R.whichCell = whichCell;

if verbose
    disp('done with everything');
    disp(['We missed ' num2str(nMiss) ' probable events across all datasets']);
    disp([' with a mean error of ' num2str(nanmean(meanErr)) ' seconds for those we got']);
    disp('congratulations');
end

