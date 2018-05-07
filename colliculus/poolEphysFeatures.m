%% params
localRoot = 'C:\DATA\Spikes\';

overwrite = true;
verbose = true;

acg_time = 1; % seconds
acg_binsize = 0.0005; % seconds
isi_time = 1; % seconds
isi_binsize = 0.001; % seconds
pspect_binsize = 0.002; % seconds
ss_fact = 10; % factor for subsampling the spectrum
% (with our choice of tapers we don't get much frequency
% resolution anyway, so might as well save space)
psparams.Fs = 1/pspect_binsize; % further parameters for power spectra
psparams.pad = 0;
psparams.fpass = [0 200];
psparams.tapers = [7 12];
psparams.trialave = true;

clear db
ephys_celltypes_db

%% run
feats = struct;
labs = {'TtP','FWHM','FW3M','Grad(0.06)','Grad(0.4)','Peak/Trough','Capacitative', ...
    'Refractory time','Mean lag','ACG Peakiness','ACG half-width','ACG 2/3-width', ...
    'CV ISI','Mode ISI','Mean rate','Spontaneous rate'};
ephysFeats = [];    ACGfeats = [];  ISIfeats = [];
% freqFeats = [];
ACG = [];           ACG_bins = [];  
distISI = [];       pxx = [];
nSpks = [];         Amp = [];       Depth = [];
allWFs = [];        whichCell = []; %isGood = [];

for k = 1:length(db)
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        if verbose
            disp(['Now: ' db(k).mouse_name ' on ' db(k).date ', ephys_' thisTag])
        end
        % load everything for that probe
        [dsetFolders, dataDir, alnDir, alfDir] = ...
            expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        if ~isempty(db(k).ksRoot)
            ksDir = [db(k).ksRoot dsetFolders '\sorting\'];
        else
            ksDir = [dataDir '\sorting\'];
        end
        
        saveFolder = [localRoot dsetFolders '\'];
        if ~exist(saveFolder,'dir')
            mkdir(saveFolder)
        elseif exist([saveFolder 'spike_feats.mat'],'file') && ~overwrite
            continue
        end
        
        aln = loadAlignments(alnDir,thisTag,db(k).tlExp,'noise',db(k).noiseExp);
        spks = loadNeuralData(ksDir,dataDir);
        
        isneu = ismember(spks.cids,unique(spks.clu(spks.isNeuron)));
        [~,max_site_temps] = max(max(abs(spks.tempsUnW),[],2),[],3); % not zero-indexed
        max_site = max_site_temps(spks.cluTemps(spks.cids+1)+1); % need to add 1 for indexing, both times
        cluSites = max_site(isneu);
        
        cluDepth = clusterAverage(spks.clu(spks.isNeuron),spks.spikeDepths(spks.isNeuron));
        cluDepth = cluDepth(:);
        cluAmps = clusterAverage(spks.clu(spks.isNeuron),spks.spikeAmps(spks.isNeuron))*spks.gain;
        
        goodCluID = spks.cids(spks.cgs==2);
        allCluID = spks.cids(isneu);
        
        cluNspk = zeros(sum(isneu),1);
        for iClu = 1:sum(isneu)
            cluNspk(iClu) = sum(spks.clu == allCluID(iClu));
        end
        
        %  optogenetics information
        if db(k).laserExp % we want to exclude laser trials from ACG and CCG analysis?
            lasOnset = []; lasOffset = []; % or not, you decide
            
            for ii = db(k).laserExp
                onsetName = sprintf('mpep_%d_onsets_in_timeline_%d.npy',ii,db(k).tlExp);
                offsetName = sprintf('mpep_%d_offsets_in_timeline_%d.npy',ii,db(k).tlExp);
                lasOnset = [lasOnset; readNPY([alnDir onsetName])];
                lasOffset = [lasOffset; readNPY([alnDir offsetName])];
            end
            laserTimes = [lasOnset lasOffset];
        end
        
        %         lasIndex = times2frames(laserTimes,tl.rawDAQTimestamps); % convert to logical indexing
        %
        %         % TO DO: check for ectopic lasers (ones not given in the alignment info)
        %         [lmin, lmax] = range(tl.raqDAQData(lasIndex,lasChannel));
        %         lasThresh = diff([lmin lmax])/2;
        
        % get spontaneous times
        beh = loadALF(alfDir,'spontaneous','sparseNoise');
        if isfield(beh,'spontaneous')
            spontTimes = beh.spontaneous.intervals(diff(beh.spontaneous.intervals,[],2) > 60,:);
            spontTimes = spontTimes - aln.tag2ref(2);
        elseif isfield(beh,'sparseNoise')
            spontTimes = [max(beh.sparseNoise.times) + 60, max(spks.st)];
            spontTimes = spontTimes - aln.tag2ref(2);
        else % assuming all recordings without sparse noise are sylvia's wheel-running expts
            spontTimes = [min(spks.st) max(spks.st)];
        end
        
        %% load or approximate the dorsal surface of SC
        bfname =[alfDir thisTag '\borders_' thisTag '.tsv'];
        snname = [saveFolder 'sparse_noise_RFs.mat'];
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
        
        %% choose which cells to keep
        % we only want cells in the mid-brain, as they're much more interesting
        mesoCells = cluDepth < scTop;% & cluDepth(:) > scBottom;
        mesoCells = mesoCells & cluNspk(:) >= 800; % got to have spikes
        
        mesoCluID = spks.cids(mesoCells);
        goodMeso = mesoCells & spks.cgs(:) == 2; % and be labelled as 'good'
        
        if verbose, disp(['  found ' num2str(sum(mesoCells)) ' cells in midbrain']); end
        if verbose, disp(['  of which ' num2str(sum(goodMeso)) ' are ''good''']); end
        if ~any(goodMeso), continue; end
        
        nClu = sum(goodMeso);
        goodMesoCIDs = spks.cids(goodMeso); % these are the cells we use
        mesoSites = cluSites(goodMeso);
        
        % remember, both cluID and cluTemps are zero-indexed!!
        cluWF = spks.tempsUnW(spks.cluTemps(goodMesoCIDs+1)+1,:,:);
        
        % first column is cluster ID, second is tag, third is db index
        % to do: find a better indexing system!
        whichCell = [whichCell; goodMesoCIDs(:) ones(nClu,1)*t ones(nClu,1)*k];
        nSpks = [nSpks; cluNspk(goodMeso)];
        Amp = [Amp; cluAmps(goodMeso)];
        Depth = [Depth; cluDepth(goodMeso), cluDepth(goodMeso) - scTop, scTop - cluDepth(goodMeso) - scBottom];
        
        %% get waveforms and features
        if verbose
            tic;
            disp('Extracting waveforms ...');
        end
        %         gwfparams.dataDir = dataDir; % for getting mean waveforms
        %         gwfparams.ksDir   = ksDir;     % takes a loooong time
        %         gwfparams.fileName = spks.dat_path;
        %         gwfparams.dataType = spks.dtype;
        %         gwfparams.nCh = spks.n_channels_dat;
        %         gwfparams.wfWin = [-31 35];
        %         gwfparams.nWf = 1000;
        %         gwfparams.spikeTimes = round(stGood*spks.sample_rate);
        %         gwfparams.spikeClusters = scGood;
        %         gwfparams.nMaxCh = [1 1];
        %         gwfparams.maxSite = goodSites;
        
        mainTemplates = zeros(nClu,size(cluWF,2));
        for iCell = 1:nClu, iSite = mesoSites(iCell);
            mainTemplates(iCell,:) = cluWF(iCell,:,iSite);
        end
        allWFs = [allWFs; mainTemplates];
        
        %         wfMainChan = getWaveFormsMA(gwfparams); % mean and first 1000 WFs on max-amplitude channel
        %         allWF = permute(wfMainChan.waveForms,[4 2 1 3]); % only needs to be 3D
        %         meanWF = permute(wfMainChan.waveFormsMean,[3 1 2]); % also, makes it easier to plot
        %         disp([' ... done at ' num2str(toc) ' seconds'])
        
        % calculate features on 'good' units
        if verbose, disp('... and computing some features ...'); end
        sfMeanWF = spikeFeatures(mainTemplates',spks.sample_rate);
        if verbose, disp([' ... done at ' num2str(toc) ' seconds']); end
        
        save([saveFolder 'spike_feats.mat'],'sfMeanWF')
        
        %         TP = [TP; sfMeanWF.TP]; FWXM = [FWXM; sfMeanWF.FWXM]; Grad = [Grad; sfMeanWF.Grad];
        %         HRatio = [HRatio; sfMeanWF.Amps]; I_cap = [I_cap; sfMeanWF.I_cap];
        
        ephysFeats = [ephysFeats; sfMeanWF.TP sfMeanWF.FWXM sfMeanWF.Grad sfMeanWF.Amps sfMeanWF.I_cap];
        
        %% ACG and CCG
        if verbose, disp('Computing autocorrelations ...'); end
        
        nbin = 2000*acg_time;
        dsetACGs = zeros(nClu,nbin + 1);
        dsetACGbins = zeros(nClu,nbin +1);
        for iCell = 1:nClu
            thesest = spks.st(spks.clu == goodMesoCIDs(iCell));
            [acg, bins] = CCG_mine(thesest,1,acg_binsize,nbin,1,[],'hz');
            bins = bins/1000;
            dsetACGs(iCell,:) = acg(bins>=0);
            dsetACGbins(iCell,:) = bins(bins>=0);
        end
        if verbose, disp([' ... done at ' num2str(toc) ' seconds']); end
        ACG = [ACG; dsetACGs];
        ACG_bins = [ACG_bins; dsetACGbins];
        
        af = ACGfeatures(dsetACGs, dsetACGbins);
        ACGfeats = [ACGfeats; af.t_ref af.t_mean af.peakiness af.PW2 af.PW3];
        
        %% ISI
        if verbose, disp('Getting ISI distributions ...'); end
        nbin = isi_time/isi_binsize;
        pIsi = nan(nClu,nbin);
        CVisi = zeros(nClu,1);
        modeISI = zeros(nClu,1);
        rMean = zeros(nClu,1);
        rMeanSpont = zeros(nClu,1);
        for iCell = 1:nClu
            thesest = spks.st(spks.clu == goodMesoCIDs(iCell));
            isi = diff(thesest);
            spontisi = isi(logical(WithinRanges(thesest(1:end-1),spontTimes)));
            isiBin = 0:isi_binsize:isi_time;
            p = histcounts(isi,isiBin);
            pIsi(iCell,:) = p;
            CVisi(iCell) = std(isi)./mean(isi);
            modeISI(iCell) = mode(isi);
            rMean(iCell) = mean(1./isi);
            rMeanSpont(iCell) = mean(1./spontisi);
        end
        distISI = [distISI; pIsi];
        if verbose, disp([' ... done at ' num2str(toc) ' seconds']); end
        ISIfeats = [ISIfeats; CVisi modeISI rMean rMeanSpont];
        
        %         %% Visual response
        %         if exist(snname,'file')
        %             load(snname)
        %         else
        %             visFeats = [visFeats; nan(nClu,1)];
        %             continue
        %         end
        %         if verbose, disp('Getting visual modulation ...'); end
        %         moddir = zeros(nClu,1);
        %         for iCell = 1:nClu
        %             dis = snrf.neur_ID == goodMesoCIDs(iCell);
        %             [~, ind1] = max(max(abs(snrf.neur_rfmap{dis}),[],2));
        %             [~, ind2] = max(max(abs(snrf.neur_rfmap{dis}),[],1));
        %             moddir(iCell) = sign(snrf.neur_rfmap{dis}(ind1,ind2))*snrf.neur_rfstats(dis).peakZscore;
        %         end
        %         visFeats = [visFeats; moddir];
        %         if verbose, disp([' ... done at ' num2str(toc) ' seconds']); end
        
        %% Rhythmicity
        if verbose, disp('Estimating power spectra ...'); end
        
        erry10_min = 0:600:(spks.nSampDat/spks.sample_rate);
        for iCell = 1:nClu
            thesest = spks.st(spks.clu == goodMesoCIDs(iCell));
            [ba,bins] = timestampsToBinned(thesest,erry10_min,pspect_binsize,[0 600]);
            [S, F, R] = mtspectrumpb(ba',psparams);
            pxx = [pxx; S(1:ss_fact:end)'];
        end
        
        if verbose, disp([' ... done at ' num2str(toc) ' seconds']); end
        
    end
end

%% package
feats.features = [ephysFeats ACGfeats ISIfeats];

feats.labels = labs;            feats.whichCell = whichCell;
feats.amplitude = Amp;          feats.nSpks = nSpks;
feats.depth = Depth;
feats.ACG = ACG;                feats.ACG_bins = ACG_bins(1,:);
feats.distISI = distISI;        feats.ISI_bins = isiBin(1:end-1);
feats.powerSpectra = pxx;       feats.F_bins = F(1:ss_fact:end);
feats.spikeShapes = allWFs;     feats.datasets = db;

if ~exist([localRoot 'SC_celltypes\all_ephys_features.mat'],'file') || overwrite
    save([localRoot 'SC_celltypes\all_ephys_features.mat'],'feats')
end

disp('done.')
