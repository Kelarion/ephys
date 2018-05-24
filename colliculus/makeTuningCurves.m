%% params
localRoot = 'C:\DATA\Spikes\';

overwrite = true;
verbose = true;

bs = 0.001; % bin size
vis_win = [-0.05 0.2]; % window for each kind of event
mov_win = [-0.2 0.2];
rwd_win = [-0.1 0.2];

clear db
ephys_bilateral_db

%%
for k = 1:length(db)
    [~,~, alnDir, alfDir, blockDir] = ...
        expDirs(db(k).mouse_name,db(k).date,db(k).tags{1},db(k).dataServer);
    beh = loadALF(alfDir,'sparseNoise','wheelMoves');
    blk = loadBlocks(blockDir,db(k).tlExp,'cw',db(k).cwExp,'pas',db(k).passiveExp);
    
    % behavioral info
    left_cont = [beh.cwStimOn.contrastLeft; beh.passiveStimOn.contrastLeft];
    right_cont = [beh.cwStimOn.contrastRight; beh.passiveStimOn.contrastRight];
    
    choice = beh.cwResponse.choice;
    feedback = beh.cwFeedback.type;
            
    wheelPos = beh.wheel.position; % wheel info
    wheelVel = beh.wheel.velocity;
    moveTypes = beh.wheelMoves.type;
    
    pas_conds = [blk.pas.trial.condition];
    pasXPos = [pas_conds.distBetweenTargets] /2; % why do I need to do this
    cwXPos = blk.cw.parameters.distBetweenTargets/2;
    cwYPos = blk.cw.parameters.targetAltitude;
    cwSigma = blk.cw.parameters.cueSigma(1);
    stimPos = [ones(length(beh.cwStimOn.contrastLeft),1)*cwXPos; pasXPos(:)];
    
    iscw = [true(length(beh.cwStimOn.times),1); false(length(beh.passiveStimOn.contrastLeft),1)];
    
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
        
        % get times for this probe
        stimTimes = -aln.tag2ref(2) + [beh.cwStimOn.times; beh.passiveStimOn.times];
        beepTimes = -aln.tag2ref(2) + [beh.cwGoCue.times];
        choiceTimes = -aln.tag2ref(2) + beh.cwResponse.times;
        feedbackTimes = -aln.tag2ref(2) + beh.cwFeedback.times;
        moveTimes = -aln.tag2ref(2) + beh.wheelMoves.intervals;
        
        % wheel times
        wt =  -aln.tag2ref(2) + linspace(beh.wheel.timestamps(1,2),beh.wheel.timestamps(2,2),length(beh.wheel.position));
        
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
        
        cwstart = aln.cw2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
        cwend = max(stimTimes(iscw));
        passtart = aln.pas2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
        pasend = max(stimTimes(~iscw));
        
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
        
        load(snname)
        %% select neurons for analysis 
        inclCells = cluDepth(:) < scTop;  % in the sSC
        inclCells = inclCells & cluNspk(:) >= 9000; % with enough spikes
%         inclCells = inclCells & snrf.neurHasRF(:); % with a visual receptive field
        
        inclCID = spks.cids(inclCells);
        nCells = length(inclCID);
        
        %% assess visual responses
        if verbose, disp(['Making plots for ' num2str(nCells) ' cells...']); end
        
        figDir = [localRoot '\' db(k).mouse_name '\' db(k).date '\visualTuning\'];
        if ~exist(figDir,'dir')
            mkdir(figDir)
        end
        
        whichVar = [right_cont left_cont];
        figure('units','normalized','position',[0.0370 0.2231 0.8786 0.6389],'visible','off');
        for iclu = 1:nCells
            fname = sprintf('Cell_%d_ephys_%s_vis.png', inclCID(iclu),thisTag);
            if ~exist([figDir fname],'file') || overwrite
                [~,bins,~,~,~,ba] = psthAndBA(spks.st(spks.clu == inclCID(iclu)),stimTimes,vis_win,bs);
                
                tuningcurve_axes
                params.whichAx = ax([1:3,5,6]);
                params.respWin = (bins>=0);
                params.varlabs = {'Right contrast','Left contrast'};
                params.majorVar = 1;
                
                plot2DTuning(ba,whichVar,bins,params)
                ax(6).XTick = [];
                xlabel(ax(6),'')
                title(ax(5),sprintf('Cell %d, %s %s %s',inclCID(iclu),db(k).mouse_name,db(k).date,thisTag))
                
                print(gcf,[figDir fname],'-dpng')
            end
        end
        close(gcf)
        
        %% other psth go here
        
        figure('units','normalized','position',[0.0370 0.2231 0.8786 0.6389],'visible','off');
        for iclu = 1:nCells
            % choice-locked, by choice made and correct choice
            whichVar = [sign(-tan(choice-3)),sign(left_cont(iscw)-right_cont(iscw))];
            figDir = [localRoot '\' db(k).mouse_name '\' db(k).date '\choiceTuning\'];
            if ~exist(figDir,'dir')
                mkdir(figDir)
            end
            fname = sprintf('Cell_%d_ephys_%s_choice.png', inclCID(iclu),thisTag);
            if ~exist([figDir fname],'file') || overwrite
                [~,bins,~,~,~,ba] = psthAndBA(spks.st(spks.clu == inclCID(iclu)),choiceTimes,mov_win,bs);
                
                tuningcurve_axes
                params.whichAx = ax([1:3,5,6]);
                params.respWin = (bins<=0);
                params.varlabs = {'Choice made','Correct choice'};
                params.majorVar = 1;
                params.useDiff = false;
                
                plot2DTuning(ba,whichVar,bins,params)
                ax(6).XTick = [];
                xlabel(ax(6),'')
                title(ax(5),sprintf('Cell %d, %s %s %s',inclCID(iclu),db(k).mouse_name,db(k).date,thisTag))
                
                print(gcf,[figDir fname],'-dpng')
            end
            
            % feedback-locked, by outcome and choice made
            whichVar = [feedback,sign(-tan(choice-3))];
            figDir = [localRoot '\' db(k).mouse_name '\' db(k).date '\feedbackTuning\'];
            if ~exist(figDir,'dir')
                mkdir(figDir)
            end
            fname = sprintf('Cell_%d_ephys_%s_feedback.png', inclCID(iclu),thisTag);
            if ~exist([figDir fname],'file') || overwrite
                [~,bins,~,~,~,ba] = psthAndBA(spks.st(spks.clu == inclCID(iclu)),feedbackTimes,rwd_win,bs);
                
                tuningcurve_axes
                params.whichAx = ax([1:3,5,6]);
                params.respWin = (bins>=0);
                params.varlabs = {'Feedback','Choice made'};
                params.majorVar = 1;
                params.useDiff = true;
                
                plot2DTuning(ba,whichVar,bins,params)
                
                ax(6).XTick = [];
                xlabel(ax(6),'')
                title(ax(5),sprintf('Cell %d, %s %s %s',inclCID(iclu),db(k).mouse_name,db(k).date,thisTag))
                
                print(gcf,[figDir fname],'-dpng')
            end
        end
        
    end
end


