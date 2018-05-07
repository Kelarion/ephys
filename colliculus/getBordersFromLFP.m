figDir = 'C:\DATA\figs\SC_project\LFPborders\';
localRoot = 'C:\DATA\Spikes\';

clear db
ephys_celltypes_db
%% 
% k = 11;
% thisTag = db(k).tags{1};
for k = 7:10
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        %% load everything
        dsetFolders = [db(k).mouse_name '\' db(k).date '\'];
        
        [ksFolders, dataDir, alnDir, blkDir, alfDir] = ... 
            expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [db(k).ksRoot ksFolders];
        
        saveFolder = [localRoot dsetFolders 'ephys_' thisTag '\'];
        
        % load timeline and blocks into a struct
        blks = loadBlocks(blkDir,db(k).tlExp,'cw',db(k).cwExp,'passive',db(k).passiveExp);
        
        % load alignment info
        aln = loadAlignments(alnDir,thisTag,db(k).tlExp, ... 
            'cw',db(k).cwExp,'noise',db(k).noiseExp,'pas',db(k).passiveExp);
        
        % get sparse noise info
        if exist([saveFolder 'sparse_noise_RFs.mat'],'file')
            tmp = load([saveFolder 'sparse_noise_RFs.mat']);
            snrf = tmp.snrf;
            try
                stimTimes = snrf.stimTimes_local;
            catch
                stimTimes = readNPY([alfDir 'sparseNoise.times.npy']) - aln.tag2ref(2);
            end
            stimPositions = snrf.stimPosition;
        elseif exist(alfDir,'dir')
            stimTimes = readNPY([alfDir 'sparseNoise.times.npy']) - aln.tag2ref(2);
            stimPositions = readNPY([alfDir 'sparseNoise.positions.npy']);
        else
            disp(['skipping ' db(k).mouse_name '_' db(k).date '_' thisTag ', no alf'])
            continue
        end
        
        % load neural data
        spks = loadNeuralData(ksDir,dataDir); 
        
        liveChans = find(ismember(0:384,spks.channel_map)); % indices of the recording channels
        
        stAll = spks.st(spks.isNeuron);
        spikeDepthsAll = spks.spikeDepths(spks.isNeuron);
        spikeAmpsAll = spks.spikeAmps(spks.isNeuron);
        
        %% LFP things
        clear bord sc scupper sclower whichsc
        bf =[alfDir thisTag '\borders_' thisTag '.tsv'];
        if exist(bf,'file')
            bord = readtable(bf ,'Delimiter','\t', 'FileType', 'text');
            sc = contains(bord.acronym,'SC');
            scupper = max(spks.ycoords) - bord.upperBorder(sc);
            sclower = max(spks.ycoords) - bord.lowerBorder(sc);
            whichsc = bord.acronym(sc);
        end
        
        % load LFP
        foo = split(spks.dat_path,'.ap');
        lfpName = [foo{1} '.lf.bin'];
        lfpnamestruct = dir([dataDir lfpName]);
        dataTypeNBytes = numel(typecast(cast(0, spks.dtype), 'uint8')); % determine number of bytes per sample
        
        nSampLFP = lfpnamestruct.bytes/(spks.n_channels_dat*dataTypeNBytes);
        lfpFs = 2500; % assumption for imec
        maxTimeForCorr = 600; % seconds; keep it < ~10 minutes or the array gets too huge
        
        lfpFile = memmapfile([dataDir lfpName], 'Format',  ...
                {spks.dtype, [spks.n_channels_dat nSampLFP], 'x'});
        
        % get start and end of each block in ephys time
        lastCWtrial = aln.cw2tl(2) + blks.cw.trial(end).trialEndedTime; % this is in 'timeline' time
        if isempty(lastCWtrial)
            lastCWtrial = aln.cw2tl(2) + blks.cw.trial(end-1).trialEndedTime + 10; 
        end
        lastPassiveTrial = aln.pas2tl(2) + blks.passive.trial(end).trialEndedTime; % in 'timeline' time
        
        tl2ephys = aln.tl2ref(2) - aln.tag2ref(2); % lag between tl time and ephys time
        
        cwTimes = tl2ephys + [aln.cw2tl(2) lastCWtrial]; % ChoiceWorld start and end in ephys time
        spontTimes = tl2ephys + [lastCWtrial aln.noise2tl(2)]; % Spontaneous start and end in ephys time
        noiseTimes = tl2ephys + [aln.noise2tl(2) aln.pas2tl(2)];
        passiveTimes = tl2ephys + [aln.pas2tl(2) lastPassiveTrial]; % in ephys time
        
        cwInds = round(cwTimes*lfpFs); % convert to LFP samples
        spontInds = round(spontTimes*lfpFs);
        noiseInds = round(noiseTimes*lfpFs);
        passiveInds = round(passiveTimes*lfpFs);
        
        blockInds = [cwInds; spontInds; noiseInds];
        blockTimes = [cwTimes; spontTimes; noiseTimes];
        blockNames = {'ChoiceWorld','Darkness','Noise'};
        nBlocks = size(blockInds,1);
        
        % load or compute LFP receptive field information
        if exist('snrf','var')
            if isfield(snrf,'stimTimes_local')
                timeCourse = snrf.lfp_timecourse;
                tcChans = snrf.computed_channels;
            else
                [timeCourse,~,tcChans] = snrfLFP(lfpFile,stimPositions,stimTimes,liveChans);
            end
        else
%             tlfp = (0:nSampLFP-1)/lfpFs;
            [timeCourse,~,tcChans] = snrfLFP(lfpFile,stimPositions,stimTimes,liveChans);
        end
        [m, ind] = max(-timeCourse,[],2);
        peakTime = ind(m == max(m));
        depthResponse = zscore(timeCourse(:,peakTime));
        respChan = tcChans(thresholdMinDur(-depthResponse, 1.5,3));
        infSC = spks.ycoords(max(respChan)); % border from visual response
        
        %% compute LFP correlations
        cormat = zeros(length(liveChans),length(liveChans),nBlocks);
        disp('Computing LFP correlations ...')
        tic
        for iblock = 1:nBlocks
            t0 = blockInds(iblock,1);
            whenCompute = t0:min([blockInds(iblock,2) t0+round(maxTimeForCorr*lfpFs)]);
            
            mySlice = double(lfpFile.Data.x(liveChans,whenCompute))'; % start at halfway
            medianSubtractedLFP = mySlice - median(mySlice')';
            cormat(:,:,iblock) = corrcoef(medianSubtractedLFP);
        end
        disp(['    ... done in ' num2str(toc) ' seconds'])
        
        rows = findDiagBlocks(cormat(:,:,3),3); % finds row of the block diagonals
        infCtx = spks.ycoords(max(rows)); % border from correlations
        %% plot
        [x, y] = meshgrid(spks.ycoords,spks.ycoords);

        fig = figure('Units','Normalized','Position',[0.0953 0.0417 0.6177 0.8796],'visible','off');
        for iAx = 1:nBlocks
            ax1 = subplot(3,3,1 + 3*(iAx-1));
            h = pcolor(x,y,cormat(:,:,iAx)); h.EdgeColor = 'none';
            set(ax1,'Xdir', 'reverse'); 
            ylabel('Depth from bottom')
            title(ax1,['Correlations during ' blockNames{iAx}])
            hold(ax1,'on'); 
            plot(xlim,[infCtx infCtx],'-.','color','k')
            plot([infCtx infCtx],ylim,'-.','color','k')
            text(2000,infCtx,'SC/Ctx (from corrs)','verticalalignment','top','color','w')
            
            ax2 = subplot(3,3,2 + 3*(iAx-1));
            stimwin  = linspace(-0.05,0.2,size(timeCourse,2));
            thesechns = spks.ycoords(ismember(spks.channel_map,tcChans));
            imagesc(ax2,stimwin,thesechns,flipud(timeCourse))
            hold(ax2,'on'); plot([0 0], ylim,'--','color','k')
            plot(xlim,[infSC infSC],'-.','color','k')
            text(0.05,infSC,'SC/Ctx (from visual)','verticalalignment','bottom')
            hold(ax2,'off');
            title(ax2,'LFP visual response')
            
            ax3 = subplot(3,3,3 + 3*(iAx-1));
            bt = [blockTimes(iAx,1), ...
                min([blockTimes(iAx,2) blockTimes(iAx,1)+maxTimeForCorr])];
            inBlock = stAll >= bt(1) & stAll <= bt(2);
            plotDriftmap(stAll(inBlock),spikeAmpsAll(inBlock),max(spks.ycoords)- spikeDepthsAll(inBlock))
            set(ax3,'Ydir','reverse')
            ylabel(''); ax3.YTick = [];
            xlabel(''); ax3.XTick = [];
            ax3.YLim = [min(spks.ycoords) max(spks.ycoords)];
            ax3.XLim = [min(stAll(inBlock)) max(stAll(inBlock))];
            hold(ax3,'on'); plot(xlim,[infCtx infCtx],'--','color','b','linewidth',2)
            axis(ax3,'off')
            title(ax3,['Spikes during ' blockNames{iAx}])
            
            if exist('scupper','var')
                j = 1;
                for ii = unique([scupper sclower])
                    hold(ax1,'on'); hold(ax2,'on'); hold(ax3,'on')
                    plot(ax1,[ii ii],ax1.YLim,'--','linewidth',1,'color','k')
                    plot(ax1,ax1.XLim,[ii ii],'--','linewidth',1,'color','k')
                    plot(ax3,ax3.XLim,[ii ii],'--','color','g','linewidth',2)
                    plot(ax2,ax2.XLim,[ii ii],'--','color','k','linewidth',2)
                    text(ax2,0.6,ii,'Known surface','verticalalignment','top')
                    hold(ax1,'off'); hold(ax2,'off'); hold(ax3,'off')
                    j = j+1;
                end
            end
        end
        fname = sprintf('%s_%s_ephys_%s_LFP_correlations.png',db(k).date,db(k).mouse_name,thisTag);
        print(fig,[figDir fname],'-djpeg')
        close(fig)
        
%         fig2 = figure('Units','normalized','position',[0.0224 0.0380 0.7589 0.8833],'visible','off');
%         j = 0;
%         for row = 1:nBlocks
%             for col = 1:nBlocks
%                 j = j+1;
%                 if col>row 
%                     continue
%                 elseif col == row
%                     ax1 = subplot(nBlocks,nBlocks,j);
%                     h = pcolor(x,y,cormat(:,:,col)); h.EdgeColor = 'none';
%                     set(ax1,'Ydir', 'reverse');
%                     title(ax1,['Correlations during ' blockNames{col}])
%                     colorbar
%                 else
%                     ax1 = subplot(nBlocks ,nBlocks,j);
%                     h = pcolor(x,y,diff(cormat(:,:,[col row]),[],3)); h.EdgeColor = 'none';
%                     set(gca,'Ydir', 'reverse');
%                     title(ax1,[blockNames{row} ' minus ' blockNames{col}])
%                     colorbar
%                 end
%                 if col==1, ylabel('Depth from surface'); end
%                 if exist('scupper','var')
%                     for ii = unique([scupper sclower])
%                         hold(ax1,'on'); 
%                         plot(ax1,[ii ii],ax1.YLim,'--','linewidth',1,'color','k')
%                         plot(ax1,ax1.XLim,[ii ii],'--','linewidth',1,'color','k')
%                     end
%                 end
%             end
%         end
%         fname = sprintf('%s_%s_ephys_%s_diff_LFPcorrs.png', ... 
%             db(k).date,db(k).mouse_name,thisTag);
%         print(fig2,[figDir fname],'-djpeg')
%         close(fig2)
        
        disp(['finished with ' db(k).mouse_name '_' db(k).date ', ' thisTag])
    end
end
%% MUA correlations
% n_depth_groups = 100; % from Andy
% spike_binning = 0.1; % seconds
% 
% inBlock = stAll > blockTimes(3,1) & stAll < blockTimes(3,2);
% blockSpks = spks.st(inBlock);
% 
% depth_group_edges = linspace(0,max(spks.ycoords),n_depth_groups+1);
% depth_group = discretize(templateDepths,depth_group_edges);
% depth_group_centers = depth_group_edges(1:end-1)+(diff(depth_group_edges)/2);
% unique_depths = 1:length(depth_group_edges)-1;
% 
% corr_edges = blockSpks(1):spike_binning:blockSpks(end);
% corr_centers = corr_edges(1:end-1) + diff(corr_edges);
% 
% binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
% for curr_depth = 1:length(unique_depths)
%     binned_spikes_depth(curr_depth,:) = histcounts(blockSpks( ...
%         ismember(spikeTempsAll(inBlock),find(depth_group == unique_depths(curr_depth)))), ...
%         corr_edges);
% end
% 
% mua_corr = corrcoef(binned_spikes_depth');
% mua_corr(all(isnan(mua_corr),2),:) = [];
% mua_corr(:,all(isnan(mua_corr),1)) = [];


