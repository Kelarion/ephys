localRoot = 'C:\DATA\Spikes\';
overwrite = 0;
verbose = 1;
lfp_thresh = 1.2; % stronger response than this counted as 

clear db
ephys_celltypes_db
%% 
for k = 1:length(db)
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        clear snrf
        % load data
        dsetFolders = [db(k).mouse_name '\' db(k).date '\'];
        [ksFolders, dataDir, alnDir, infoDir, alfDir] = ... 
            expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [db(k).ksRoot ksFolders];
        
        saveFolder = [localRoot dsetFolders 'ephys_' thisTag '\'];
        if exist([saveFolder 'sparse_noise_RFs.mat'],'file') && ~overwrite
            continue
        elseif verbose
            disp(['Now: ' db(k).mouse_name ' on ' db(k).date ', ephys_' thisTag])
        end
        
        aln = loadAlignments(alnDir,thisTag,db(k).tlExp, 'noise',db(k).noiseExp);
        
        if exist([alfDir 'sparseNoise.times.npy'],'file')
            snrf.stimTimes_local = readNPY([alfDir 'sparseNoise.times.npy']) - aln.tag2ref(2);
            snrf.stimPosition = readNPY([alfDir 'sparseNoise.positions.npy']);
        else
            pname = sprintf('%s_%d_%s_parameters.mat',db(k).date,db(k).noiseExp,db(k).mouse_name);
            hwname = sprintf('%s_%d_%s_hardwareInfo.mat',db(k).date,db(k).noiseExp,db(k).mouse_name);
            try
                pmtr = loadVar([infoDir num2str(db(k).noiseExp) '\' pname],'parameters');
                hw = loadVar([infoDir num2str(db(k).noiseExp) '\' hwname],'myScreenInfo');
                [st, sp] = mpep2stimTimes(aln.noise2tl(2),pmtr.Protocol,hw,aln.tl2ref);
                snrf.stimTimes_local = st{1} - aln.tag2ref(2); 
                snrf.stimPosition = sp{1}; 
            catch
                continue
            end
        end
        
        spks = loadNeuralData(ksDir,dataDir); % loads all neural data
        liveChans = find(ismember(0:384,spks.channel_map)); % indices of the recording channels
        nCells = length(unique(spks.clu(spks.isNeuron)));
        
        %% get RFs from visual noise
        % for LFP
        if verbose, disp('Getting LFP RFs...'); end
        lfpFile = memmapfile([dataDir spks.lfp_path], 'Format',  ...
                {spks.dtype, [spks.n_channels_dat spks.nSampLFP], 'x'});
        [timeCourse, lfpRF, chn] = snrfLFP(lfpFile,snrf.stimPosition,snrf.stimTimes_local,liveChans);
        
        [m, ind] = max(-timeCourse,[],2);
        peakTime = ind(m == max(m));
        depthResponse = zscore(timeCourse(:,peakTime));
        respChan = chn(thresholdMinDur(-depthResponse, lfp_thresh,3));
        
        snrf.lfp_rfmap = lfpRF;   snrf.lfp_timecourse = timeCourse;  
        snrf.responsive_channels = respChan; snrf.computed_channels = chn;
        
        % for single cells
        if verbose, disp('Getting neural RFs...'); end
        clear rfstats rfmap
        rfmap = cell(nCells,1);
        hasRF = zeros(nCells,1);
        for iCell = 1:nCells
            params.useSVD = true; params.makePlots = false;
            params.fit2Dgauss = true;
            theseSpikes = spks.st(spks.isNeuron & spks.clu == spks.cids(iCell));
            [rfmap{iCell}, rfstats(iCell), ~] = sparseNoiseRF(theseSpikes, ... 
                snrf.stimTimes_local,snrf.stimPosition,params);
            if rfstats(iCell).fit2D(1) > 0.2
                hasRF(iCell) = true;
            end
        end
        snrf.neur_rfmap = rfmap;   snrf.neur_rfstats = rfstats;   snrf.neurHasRF = hasRF;
        
        save([saveFolder 'sparse_noise_RFs.mat'],'snrf')
        
        if verbose, disp('done.'); end
        
    end
end















