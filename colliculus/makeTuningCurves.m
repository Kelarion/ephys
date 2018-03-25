%% pool everything
localRoot = 'C:\DATA\Spikes\';

overwrite = 1;
verbose = 1;

clear db
ephys_bilateral_db

%%
for k = 1:length(db)
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        
        thisExp =  [db(k).mouse_name '_' db(k).date '_' thisTag];
        if verbose, disp(['Now: ' thisExp]); end
        
        % load everything for that probe
        [dsetFolders, dataDir, alnDir, blockDir, alfDir] = ... 
            expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [dataDir '\sorting\'];
        
        saveFolder = [localRoot dsetFolders '\'];
        if ~exist(saveFolder,'dir')
            mkdir(saveFolder)
        end
        snname = [saveFolder 'sparse_noise_RFs.mat'];
        if exist(snname,'file')
            snrf = loadVar(snname,'snrf');
        else
            disp([thisExp ' doesn''t have a receptive field file, moving on'])
%             continue
        end
        
        aln = loadAlignments(alnDir,thisTag,db(k).tlExp, ... 
            'noise',db(k).noiseExp,'cw',db(k).cwExp,'pas',db(k).passiveExp);
        spks = loadNeuralData(ksDir,dataDir);
        beh = loadALF(alfDir);
        blk = loadBlocks(blockDir,db(k).tlExp,'cw',db(k).cwExp,'pas',db(k).passiveExp);
        
        % behavioral info
        left_cont = [beh.cwStimOn.contrastLeft; beh.passiveStimOn.contrastLeft];
        right_cont = [beh.cwStimOn.contrastRight; beh.passiveStimOn.contrastRight];
        beepTimes = -aln.tag2ref(2) + [beh.cwGoCue.times];
        choice = beh.cwResponse.choice;
        choiceTimes = -aln.tag2ref(2) + beh.cwResponse.times;
        stimTimes = -aln.tag2ref(2) + [beh.cwStimOn.times; beh.passiveStimOn.times];
        iscw = [true(length(beh.cwStimOn.times),1); false(length(beh.passiveStimOn.contrastLeft),1)];
        
        wheelPos = beh.wheel.position; % wheel info
        wheelVel = beh.wheel.velocity;
        moveTimes = beh.wheelMoves.intervals;
        moveTypes = beh.wheelMoves.type;
        wt =  -aln.tag2ref(2) + linspace(beh.wheel.timestamps(1,2),beh.wheel.timestamps(2,2),length(beh.wheel.position));
        
        pas_conds = [blk.pas.trial.condition];
        pasXPos = [pas_conds.distBetweenTargets] /2; % why do I need to do this
        cwXPos = blk.cw.parameters.distBetweenTargets/2;
        cwYPos = blk.cw.parameters.targetAltitude;
        cwSigma = blk.cw.parameters.cueSigma(1);
        stimPos = [ones(length(beh.cwStimOn.contrastLeft),1)*cwXPos; pasXPos(:)];
        
        % further ephys info
        isneu = ismember(spks.cids,unique(spks.clu(spks.isNeuron)));
        [~,max_site_temps] = max(max(abs(spks.tempsUnW),[],2),[],3); % not zero-indexed
        max_site = max_site_temps(spks.cluTemps(spks.cids+1)+1); % need to add 1 for indexing, both times
        cluSites = max_site(isneu);
        allCluID = spks.cids(isneu);
        
        cluNspk = zeros(sum(isneu),1);
        for iClu = 1:sum(isneu)
            cluNspk(iClu) = sum(spks.clu == allCluID(iClu));
        end
        
        cluDepth = clusterAverage(spks.clu(spks.isNeuron),spks.spikeDepths(spks.isNeuron));
        cluAmps = clusterAverage(spks.clu(spks.isNeuron),spks.spikeAmps(spks.isNeuron))*spks.gain;
        
        %% load or approximate the dorsal surface of SC
        bfname =[alfDir thisTag '\borders_' thisTag '.tsv'];
        if exist(bfname,'file')
            if verbose, disp('Reading histology data...'); end
            bord = readtable(bfname ,'Delimiter','\t', 'FileType', 'text');
            sc = contains(bord.acronym,'SC');
            scupper = bord.upperBorder(sc);
            sclower = bord.lowerBorder(sc);
            scTop = max(scupper);
            scBottom = min(sclower);
        elseif exist(snname,'file')
            if verbose, disp('Estimating SC surface...'); end
            ctxChan = wheresCortex([dataDir spks.lfp_path],'snrf',snname);
            if isempty(ctxChan)
                t0 = aln.noise2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
                if verbose, disp('No visual response, using LFP correlations...'); end
                ctxChan = wheresCortex([dataDir spks.lfp_path],'corrs',[t0 t0+400]);
                if verbose, disp(['   ... I think it''s at channel: ' num2str(ctxChan)]); end
            end
            scTop = spks.ycoords(ctxChan);
            scBottom = 0;
        elseif exist('aln','var')
            t0 = aln.noise2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
            if verbose, disp('Estimating SC surface through LFP correlations...'); end
            ctxChan = wheresCortex([dataDir spks.lfp_path],'corrs',[t0 t0+400]);
            if verbose, disp(['   ... I think it''s at channel: ' num2str(ctxChan)]); end
            scTop = spks.ycoords(ctxChan);
            scBottom = 0;
        else % if this is the case, what business does this dataset have in the db?
            % probably not much -- take it out
            scTop = max(spks.ycoords);
            scBottom = 0;
        end
        
        load(snname)
        %% select neurons for analysis 
        inclCells = cluDepth(:) < scTop & cluDepth(:) > scBottom; % in the midbrain
        inclCells = inclCells & cluNspk(:) >= 9000; % with enough spikes
        inclCells = inclCells & snrf.neurHasRF(:); % with a visual receptive field
        
        inclCID = spks.cids(inclCells);
        nCells = length(inclCID);
        
        %% assess visual responses
        
        
        
    end
end


