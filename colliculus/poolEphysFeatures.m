%% pool everything
localRoot = 'C:\DATA\Spikes\';

overwrite = 1;
verbose = 1;

acg_time = 1; % second
acg_binsize = 0.0005; % seconds
isi_time = 1; % seconds
isi_binsize = 0.001; % seconds

clear db
ephys_celltypes_db

feats = struct;
labs = {'TtP','FWHM','FW3M','Grad(0.06)','Grad(0.4)','Peak/Trough','I_cap', ...
        'Refractory time','Mean lag','ACG Peakiness','Mean firing rate', ... 
        'ACG half-width','ACG 2/3-width', ...
        'ISI CV','Mode ISI'};
ephysFeats = [];    ACGfeats = [];      ISIfeats = [];
ACG = [];           ACG_bins = [];      distISI = [];   
nSpks = [];         Amp = [];           allWFs = [];
whichCell = [];     isGood = [];            

%% 
for k = 1:length(db)
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        if verbose
            disp(['Now: ' db(k).mouse_name ' on ' db(k).date ', ephys_' thisTag])
        end
        % load everything for that probe
        dsetFolders = [db(k).mouse_name '\' db(k).date '\'];
        
        [ksFolders, dataDir, alnDir, ~, alfDir] = ... 
            expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [db(k).ksRoot ksFolders];
        
        saveFolder = [localRoot dsetFolders 'ephys_' thisTag '\'];
        if ~exist(saveFolder,'dir')
            mkdir(saveFolder)
        elseif exist([saveFolder 'spike_feats.mat'],'file') && ~overwrite
            continue
        end
        if exist(alnDir,'dir') && isfield(db,'noiseExp')
            aln = loadAlignments(alnDir,thisTag,db(k).tlExp,'noise',db(k).noiseExp);
        end
%         % load timeline for optogenetics information
%         tl = load([infoRoot expFolders num2str(db(k).tlExp) '\' tlInfo '_Timeline.mat']);
%         lasChannel = tl.inputs(strcmp({tl.hw.inputs.name},'blueLEDmonitor')).arrayColumn;
%         
%         if db(k).laserExp % we want to exclude laser trials from ACG and CCG analysis
%             lasOnset = []; lasOffset = []; % or not, you decide
%             
%             for ii = 1:length(db(k).laserExp)
%                 onsetName = sprintf('mpep_%d_onsets_in_timeline_%d.npy',ii,db(k).tlExp);
%                 offsetName = sprintf('mpep_%d_offsets_in_timeline_%d.npy',ii,db(k).tlExp);
%                 lasOnset = [lasOnset; readNPY([alignFolder onsetName])];
%                 lasOffset = [lasOffset; readNPY([alignFolder offsetName])];
%             end
%             laserTimes = [lasOnset lasOffset];
%         end
%         
%         lasIndex = times2frames(laserTimes,tl.rawDAQTimestamps); % convert to logical indexing
%         
%         % TO DO: check for ectopic lasers (ones not given in the alignment info)
%         [lmin, lmax] = range(tl.raqDAQData(lasIndex,lasChannel));
%         lasThresh = diff([lmin lmax])/2;

        %% load neural data
        spks = loadNeuralData(ksDir,dataDir);

        [~,max_site_temps] = max(max(abs(spks.tempsUnW),[],2),[],3); % not zero-indexed
        max_site = max_site_temps(spks.cluTemps(spks.cids+1)+1); % need to add 1 for indexing, both times
        cluSites = max_site(spks.cgs ~= 3);
        
        cluDepth = clusterAverage(spks.clu(spks.isNeuron),spks.spikeDepths(spks.isNeuron));
        cluAmps = clusterAverage(spks.clu(spks.isNeuron),spks.spikeAmps(spks.isNeuron))*spks.gain;
        
        goodCluID = spks.cids(spks.cgs==2);
        allCluID = spks.cids(spks.cgs ~= 0);
        
        cluNspk = zeros(sum(spks.cgs ~= 0),1);
        for iClu = 1:sum(spks.cgs ~= 0)
            cluNspk(iClu) = sum(spks.clu == allCluID(iClu));
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
        
        %% choose which cells to keep
        % we only want cells in the mid-brain, they're much more interesting anyway
        mesoCells = cluDepth < scTop & cluDepth > scBottom;
        mesoCells = mesoCells & cluNspk >= 800; % got to have spikes
        
        nMesoUnits = sum(mesoCells);
        mesoCluID = spks.cids(mesoCells);
        mesoSites = cluSites(mesoCells);
        nGoodMeso = sum(ismember(mesoCluID,goodCluID));
        
        if verbose, disp(['  found ' num2str(nMesoUnits) ' cells in midbrain']); end
        if verbose, disp(['  of which ' num2str(nGoodMeso) ' are ''good''']); end
        
        % remember, both cluID and cluTemps are zero-indexed!!
        cluWF = spks.tempsUnW(spks.cluTemps(mesoCluID+1)+1,:,:);
        
        % add to structure
        spks.mesoCluID = mesoCluID;
        spks.scTop = scTop;
        spks.scBottom = scBottom;
        % first column is cluster ID, second is dataset index, third is tag
        % to do: find a better way to index datasets!
        whichCell = [whichCell; mesoCluID' ones(nMesoUnits,1)*k ones(nMesoUnits,1)*t];
        isGood = [isGood; ismember(mesoCluID,goodCluID)'];
        nSpks = [nSpks; cluNspk(mesoCells)];
        Amp = [Amp; cluAmps(mesoCells)];
        allSpks(k+(t-1)) = spks;
        
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

        mainTemplates = zeros(nMesoUnits,size(cluWF,2));
        for iCell = 1:nMesoUnits, iSite = mesoSites(iCell);
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
        dsetACGs = zeros(nMesoUnits,nbin + 1);
        dsetACGbins = zeros(nMesoUnits,nbin +1);
        for iCell = 1:nMesoUnits
            thesest = spks.st(spks.clu == mesoCluID(iCell));
            [acg, bins] = CCG_mine(thesest,1,acg_binsize,nbin,1);
            bins = bins/1000;
            dsetACGs(iCell,:) = acg(bins>=0);
            dsetACGbins(iCell,:) = bins(bins>=0);
        end
        if verbose, disp([' ... done at ' num2str(toc) ' seconds']); end
        ACG = [ACG; dsetACGs];
        ACG_bins = [ACG_bins; dsetACGbins];
        
        af = ACGfeatures(dsetACGs, dsetACGbins);
        ACGfeats = [ACGfeats; af.t_ref af.t_mean af.peakiness af.r_mean af.PW2 af.PW3];
        
        %% ISI distribution
        if verbose, disp('Getting ISI distributions ...'); end
        nbin = isi_time/isi_binsize;
        pIsi = nan(nMesoUnits,nbin);
        CVisi = zeros(nMesoUnits,1);
        isiMax = zeros(nMesoUnits,1);
        for iCell = 1:nMesoUnits
            thesest = spks.st(spks.clu == mesoCluID(iCell));
            isi = diff(thesest);
            isiBin = 0:isi_binsize:isi_time;
            p = histcounts(isi,isiBin);
            pIsi(iCell,:) = p;
            CVisi(iClu) = mean(isi)/var(isi);
            [m,ind] = max(p);
            isiMax(iClu) = isiBin(ind);
%             p_decay = p(ind:end);
%             t_tau = isiBin(find(p_decay<m*exp(-1),1,'first')) + isiBin(ind-1);
        end
        distISI = [distISI; pIsi];
        if verbose, disp([' ... done at ' num2str(toc) ' seconds']); end
        ISIfeats = [ISIfeats; CVisi isiMax];

    end
end

% package
feats.features = [ephysFeats ACGfeats ISIfeats];       
feats.featureLabels = labs;
feats.amplitude = Amp;          feats.nSpks = nSpks;
feats.ACG = ACG;                feats.ACG_bins = ACG_bins;             
feats.distISI = distISI;        feats.ISI_bins = isiBin(1:end-1);
feats.spikeShapes = allWFs;     feats.sp = allSpks;
feats.whichCell = whichCell;    feats.isGood = logical(isGood);
feats.datasets = db;

if ~exist([localRoot 'SC_celltypes\all_ephys_features.mat'],'file') || overwrite
    save([localRoot 'SC_celltypes\all_ephys_features.mat'],'feats')
end

disp('done.')
