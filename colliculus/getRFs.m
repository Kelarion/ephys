localRoot = 'C:\DATA\Spikes\';

overwrite = false;
verbose = true;

lfp_thresh = 1.2; % threshold response when detecting sSC surface
z_thresh = 9; % threshold on peak RF Z-score
params.useSVD = true; 
params.makePlots = false;
params.fit2Dgauss = false; % we're doing this later anyway

clear db % first half of script is sensitive to structure of 'db'
ephys_bilateral_db 
%% 
% supply a 'paragon' structure, against which any existing structure is
% compared. Make sure that this file is exactly how you want all others to
% look, since extra fields will be deleted.
% paragon = [localRoot 'SS072\2016-12-02\ephys_SC\sparse_noise_RFs.mat'];
for k = 1:length(db)
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        clear snrf
        % load data
        [dsetFolders,dataDir, alnDir, alfDir,blockDir] = ... 
            expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        
        if isfield(db,'ksRoot')
            if ~isempty(db(k).ksRoot)
                ksDir = [db(k).ksRoot dsetFolders '\sorting\'];
            else
                ksDir = [dataDir '\sorting\'];
            end
        else
            ksDir = [dataDir '\sorting\'];
        end
        
        saveFolder = [localRoot dsetFolders];
        if ~overwrite && exist([saveFolder 'sparse_noise_RFs.mat'],'file')
            continue
        end
            
        spks = loadNeuralData(ksDir,dataDir); % loads all neural data
        
        aln = loadAlignments(alnDir,thisTag,db(k).tlExp, 'noise',db(k).noiseExp);
        
        if exist([alfDir 'sparseNoise.times.npy'],'file')
            stimTimes_local = readNPY([alfDir 'sparseNoise.times.npy']) - aln.tag2ref(2);
            stimPosition = readNPY([alfDir 'sparseNoise.positions.npy']);
        else
            pname = sprintf('%s_%d_%s_parameters.mat',db(k).date,db(k).noiseExp,db(k).mouse_name);
            hwname = sprintf('%s_%d_%s_hardwareInfo.mat',db(k).date,db(k).noiseExp,db(k).mouse_name);
            try
                pmtr = loadVar([blockDir num2str(db(k).noiseExp) '\' pname],'parameters');
                hw = loadVar([blockDir num2str(db(k).noiseExp) '\' hwname],'myScreenInfo');
                [st, sp] = mpep2stimTimes(aln.noise2tl(2),pmtr.Protocol,hw,aln.tl2ref);
                stimTimes_local = st{1} - aln.tag2ref(2); 
                stimPosition = sp{1}; 
            catch
                continue
            end
        end
        
        if exist([saveFolder 'sparse_noise_RFs.mat'],'file') && exist('paragon','var') 
            if verbose % save time by only fixing what is missing
                disp(['Fixing: ' db(k).mouse_name ' on ' db(k).date ', ephys_' thisTag]); 
            end
            g = load(paragon);
            snrf = loadVar([saveFolder 'sparse_noise_RFs.mat'],'snrf');
            mustHave = fieldnames(g.snrf);
            doesHave = fieldnames(snrf);
            missing = ~contains(mustHave,doesHave);
            extra = ~contains(doesHave,mustHave);
            if any(contains(mustHave(missing),'neur_ID'))
                snrf.neur_ID = spks.cids;
            end
            if any(contains(mustHave(missing),'XPos'))
                snrf.XPos = unique(stimPosition(:,2));
                snrf.YPos = -unique(stimPosition(:,1));
            end
            snrf = rmfield(snrf,doesHave(extra));
            save([saveFolder 'sparse_noise_RFs.mat'],'snrf')
            continue
        end
        if verbose
            disp(['Now: ' db(k).mouse_name ' on ' db(k).date ', ephys_' thisTag])
        end
        
        liveChans = find(ismember(0:384,spks.channel_map)); % indices of the recording channels
        nCells = length(unique(spks.clu(spks.isNeuron)));
        
        snrf.stimTimes_local = stimTimes_local;
        snrf.stimPosition = stimPosition;
        
        %% get RFs from visual noise
        % for LFP
        if verbose, disp('Getting LFP RFs...'); end
        lfpFile = memmapfile([dataDir spks.lfp_path], 'Format',  ...
                {spks.dtype, [spks.n_channels_dat spks.nSampLFP], 'x'});
        [timeCourse, lfpRF, chn] = snrfLFP(lfpFile,snrf.stimPosition,snrf.stimTimes_local,liveChans);
        
        [m, ind] = max(-timeCourse,[],2);
        peakTime = ind(m == max(m));
        depthResponse = zscore(timeCourse(:,peakTime));
        respChan = chn(thresholdMinDur(-depthResponse, lfp_thresh,2));
        
        snrf.lfp_rfmap = lfpRF;   snrf.lfp_timecourse = timeCourse;  
        snrf.responsive_channels = respChan; snrf.computed_channels = chn;
        
        % for single cells
        if verbose, disp('Getting neural RFs...'); end
        clear rfstats rfmap
        rfmap = cell(nCells,1);
        hasRF = zeros(nCells,1);
        for iCell = 1:nCells
            theseSpikes = spks.st(spks.isNeuron & spks.clu == spks.cids(iCell));
            [rfmap{iCell}, rfstats(iCell), ~] = sparseNoiseRF(theseSpikes, ... 
                snrf.stimTimes_local,snrf.stimPosition,params);
            if rfstats(iCell).peakZscore > z_thresh
                hasRF(iCell) = true;
            end
        end
        snrf.neur_rfmap = rfmap;   snrf.neur_rfstats = rfstats;   snrf.neurHasRF = hasRF;
        snrf.XPos = unique(snrf.stimPosition(:,2)); 
        snrf.YPos = -unique(snrf.stimPosition(:,1)); % it's reversed for some reason
        snrf.neur_ID = spks.cids;
        
        save([saveFolder 'sparse_noise_RFs.mat'],'snrf')
        
        if verbose, disp('done.'); end
        
    end
end





