localRoot = 'C:\DATA\Spikes\';

overwrite = true;
verbose = true;

bs = 0.01; % bin size for psth
vis_win = [-0.05 0.2]; % window for each kind of event

ccg_bs = 0.01;
ccg_win = 0.1; % +/- from 0, in sec
min_spont = 60; % minimum interval to be considered 'spontaneous', in sec

clear db
ephys_bilateral_db

% This part is strange:
% Since we are pooling across probes, CIDs will overlap between tags,
% we need to change them such that they are unique over the whole dataset.
% Ideally, we should also get back out from the transformed CIDs which
% probe and even which dataset they came from. The mapping below is what I 
% settled on, it goes from [CID, tag, dset] -> uniqueCID. I'm sure there's a
% better and simpler way, but I couldn't think of it in the moment. It's in 
% essence a conversion of base; for some reason, if the base is one greater
% than the maximum value of dset, it guarantees this mapping is invertible.
% (see testhash.m). Maybe we can just use it on [CID,tag]?
b = length(db) + 1;
wc2id = @(WC) WC*[b^2 b^1 b^0]'; % from [CID,tag,dset] -> unique CID
id2wc = @(id) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)]; % inverse

%%
TotCorrs = cell(1,4);
SigCorrs = cell(1,4);
XCorrs = cell(1,4);
for k = 1:length(db)
    [~,~, alnDir, alfDir, blockDir] = ...
        expDirs(db(k).mouse_name,db(k).date,db(k).tags{1},db(k).dataServer);
    beh = loadALF(alfDir,'spontaneous','sparseNoise');
    blk = loadBlocks(blockDir,db(k).tlExp,'cw',db(k).cwExp,'pas',db(k).passiveExp);
    
    % behavioral info
    left_cont = [beh.cwStimOn.contrastLeft; beh.passiveStimOn.contrastLeft];
    right_cont = [beh.cwStimOn.contrastRight; beh.passiveStimOn.contrastRight];
    
    choice = beh.cwResponse.choice;
    feedback = beh.cwFeedback.type;
    
%     wheelPos = beh.wheel.position; % wheel info
%     wheelVel = beh.wheel.velocity;
%     moveTypes = beh.wheelMoves.type;
    
    pas_conds = [blk.pas.trial.condition];
    pasXPos = [pas_conds.distBetweenTargets] /2; % why do I need to do this
    cwXPos = blk.cw.parameters.distBetweenTargets/2;
    cwYPos = blk.cw.parameters.targetAltitude;
    cwSigma = blk.cw.parameters.cueSigma(1);
    stimPos = [ones(length(beh.cwStimOn.contrastLeft),1)*cwXPos; pasXPos(:)];

    spontTimes = beh.spontaneous.intervals(diff(beh.spontaneous.intervals,[],2) > min_spont,:);
    stimTimes = [beh.cwStimOn.times; beh.passiveStimOn.times];
    %     spontTimes = -aln.tag2ref(2) + spont;
    %     beepTimes = -aln.tag2ref(2) + [beh.cwGoCue.times];
    %     choiceTimes = -aln.tag2ref(2) + beh.cwResponse.times;
    %     feedbackTimes = -aln.tag2ref(2) + beh.cwFeedback.times;
    %     moveTimes = -aln.tag2ref(2) + beh.wheelMoves.intervals;
    %     linspace(beh.wheel.timestamps(1,2),beh.wheel.timestamps(2,2),length(beh.wheel.position));
    
    iscw = [true(length(beh.cwStimOn.times),1); false(length(beh.passiveStimOn.contrastLeft),1)];
    
    %% pool spikes from all probes, save only relevant info
    cidAndTag = []; nSpks = [];	Amp = [];   Depth = [];
    st_all = [];    clu_all = [];
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        
        thisExp =  [db(k).mouse_name '_' db(k).date '_' thisTag];
        if verbose, disp(['Now: ' thisExp]); end
        
        % load everything for that probe
        [dsetFolders, dataDir] = expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [dataDir '\sorting\'];
        
        snrfFolder = [localRoot dsetFolders '\'];
        snname = [snrfFolder 'sparse_noise_RFs.mat'];
        if exist(snname,'file')
            snrf = loadVar(snname,'snrf');
        else
            disp([thisExp ' doesn''t have a receptive field file, moving on'])
            %             continue
        end
        
        aln = loadAlignments(alnDir,thisTag,db(k).tlExp, ...
            'noise',db(k).noiseExp,'cw',db(k).cwExp,'pas',db(k).passiveExp);
        spks = loadNeuralData(ksDir,dataDir);
                
        % further ephys info
        isneu = ismember(spks.cids,unique(spks.clu(spks.isNeuron)));
        [~,max_site_temps] = max(max(abs(spks.tempsUnW),[],2),[],3); % not zero-indexed
        max_site = max_site_temps(spks.cluTemps(spks.cids+1)+1); % need to add 1 for indexing, twice
        cluSites = max_site(isneu);
        allCluID = spks.cids(isneu);
        
        cluNspk = zeros(sum(isneu),1);
        for iClu = 1:sum(isneu)
            cluNspk(iClu) = sum(spks.clu == allCluID(iClu));
        end
        
        cluDepth = clusterAverage(spks.clu(spks.isNeuron),spks.spikeDepths(spks.isNeuron));
        cluAmps = clusterAverage(spks.clu(spks.isNeuron),spks.spikeAmps(spks.isNeuron))*spks.gain;
        
%         cwstart = aln.cw2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
%         cwend = max(stimTimes(iscw));
%         passtart = aln.pas2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
%         pasend = max(stimTimes(~iscw));
        
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
            if verbose, disp('Finding SC surface...'); end
            [ctxChan, pagchan] = wheresCortex([dataDir spks.lfp_path],'snrf',snname);
            if isempty(ctxChan)
                t0 = aln.noise2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
                ctxChan = wheresCortex([dataDir spks.lfp_path],'corrs',[t0 t0+400]);
                scTop = spks.ycoords(ctxChan);
                scBottom = 0;
            else
                scTop = spks.ycoords(ctxChan);
                scBottom = spks.ycoords(pagchan);
            end
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
        
        %% select neurons for analysis
        inclCells = cluDepth(:) < scTop;% & cluDepth(:) > scBottom;  % in the midbrain
        inclCells = inclCells & cluNspk(:) >= 9000; % with enough spikes
        %         inclCells = inclCells & snrf.neurHasRF(:); % with a visual receptive field
        
        inclCID = spks.cids(inclCells);
        nCells = length(inclCID);
       
        % Depth is [absoluteDepth, depthFromOpticLayer]
        Depth = [Depth; cluDepth(inclCells), cluDepth(inclCells) - scBottom]; 
        cidAndTag = [cidAndTag; inclCID(:) ones(nCells,1)*t ones(nCells,1)*k];
        nSpks = [nSpks; cluNspk(inclCells)];
        Amp = [Amp; cluAmps(inclCells)];
        
        % apply the bizare mapping to clu_all; I don't do it to the CIDs yet 
        % because I want to sort by tag later
        clus = [double(spks.clu(ismember(spks.clu,inclCID))) ...
            ones(sum(cluNspk(inclCells)),1)*t ones(sum(cluNspk(inclCells)),1)*k];
        
        st_all = [st_all; applyCorrection(spks.st(ismember(spks.clu,inclCID)),aln.tag2ref)];
        clu_all = [clu_all; wc2id(clus)];

    end
    %% Correlations
    whichStim = right_cont-left_cont; % a 'stimulus' is defined as bilateral visual scene
    stimTypes = unique(whichStim);
    
    nCellTot = size(cidAndTag,1);
    vis_nbin = diff(vis_win)/bs;
    ccg_nbin = round(ccg_win/ccg_bs);
    
    % define neuron order and convert to unique IDs
    [~,cluOrder] = sortrows([cidAndTag(:,2) -Depth(:,2)]); % cells organised by depth
    cid_all = wc2id(cidAndTag(cluOrder,:,:));
    
    % comodulation during stimulus presentation
    visClu = cid_all(Depth(cluOrder,2) > 0); % only look at sSC cells
    nCellVis = length(visClu);
    
    P = cell(length(unique(whichStim)),1);
    for iStim = 1:length(stimTypes) % a different permutation of stim repeat for each cell
        inds = find(whichStim == stimTypes(iStim));
        P{iStim} = nkperms(inds,nCellVis);
    end
    
    stimRespArray = zeros(nCellVis,length(stimTimes),vis_nbin);
    for iclu = 1:nCellVis
        [ba,bins] = timestampsToBinned(st_all(clu_all == cid_all(iclu,1)),stimTimes,bs,vis_win);
        stimRespArray(iclu,:,:) = ba;
    end
    
    totcorr = zeros(nCellVis,nCellVis,vis_nbin);
%     totcov = zeros(nCellVis,nCellVis,vis_nbin);
    for tau = 1:vis_nbin
        totcorr(:,:,tau) = corrcoef(stimRespArray(:,:,tau)');
%         totcov(:,:,tau) = cov(stimRespArray(:,:,tau)');
    end
    
    sigcorr = zeros(nCellVis,nCellVis,vis_nbin);
    shuffled = [];
    for iStim = 1:length(stimTypes) % shuffle array by repeats of each stimulus
        shuf = chopnotslop(stimRespArray,P{iStim},[2 1]);
        shuffled = cat(2,shuffled,shuf); 
    end
    for tau = 1:vis_nbin
        sigcorr(:,:,tau) = corrcoef(shuffled(:,:,tau)');
    end
    
    % spontaneous correlations
    [allccg,bins_ccg] = CCG(st_all,clu_all, ccg_bs, ccg_nbin,1,cid_all,'scale',spontTimes);
    
    % package
    TotCorrs{k} = totcorr;
    SigCorrs{k} = sigcorr;
    XCorrs{k} = allccg;
    
end


