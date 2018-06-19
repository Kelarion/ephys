%% Dataset
db_ephys_opticTract

%% Define folders
% protocolFolder = '\\ZSERVER.cortexlab.net\Data\trodes';
protocolFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
hardwareFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
timelineFolder = '\\ZSERVER.cortexlab.net\Data\expInfo';
subjectsFolder = '\\ZSERVER.cortexlab.net\Data\Subjects';
plotFolder = 'C:\STORAGE\OneDrive - University College London\Lab\RESULTS\OpticTract';

%% Parameters
runningSigma = 0.25; % sec, for smoothing of running speed
spikingSigma = 0.5; % sec, for smoothing spike rate

% for RFs
params.makePlots = false;
params.useSVD = true;
% params.useSVD = false;
params.countWindow = [-0.05 0.18];
params.binSize = 0.01;

RFtypes = {'Absolute'}; %, 'White', 'Black', 'Linear'};
grad = linspace(0,1,40)';
reds = [ones(40,1),grad,grad];
blues = [grad,grad,ones(40,1)];
cm = [blues; flip(reds(1:end-1,:),1)];
tmp = lines(3);
colors = [tmp(1,:); tmp(3,:); 0 0 0];

nPseudo = 20;

for k = 4:length(db)
    alignDir = fullfile(subjectsFolder, db(k).subject, db(k).date, 'alignments');
    %% Load spike data
    sp = loadAllKsDir(db(k).subject, db(k).date);
    % select sorted neurons
    % NOTE! if several electrodes were used for recording, times need to be
    % aligned! (see scripts\electrophys\extractData
    goodUnits = sp.cids(sp.cgs == 2);
    templates = findTempForEachClu(sp.clu, sp.spikeTemplates);
    
    tags = getEphysTags(db(k).subject, db(k).date);
    [~, ~, ~, ~, ~, tl, hasTimeline] = ...
        dat.whichExpNums(db(k).subject, db(k).date);
    TLexp = find(hasTimeline);
    tl = tl{1};
    
    %% Load and prepare data for receptive fields
    % data = load(fullfile(protocolFolder, subject, date, num2str(exp), ...
    %     'Protocol.mat'));
    data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_parameters.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    pars = data.parameters.Protocol;
    stimFile = str2func(strtok(pars.xfile, '.'));
    % load myScreenInfo
    load(fullfile(hardwareFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expNoise), sprintf('%s_%d_%s_hardwareInfo.mat', ...
        db(k).date, db(k).expNoise, db(k).subject)));
    myScreenInfo.windowPtr = NaN;
    
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp.name)));
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expNoise, TLexp)));
    
    noiseOn = applyCorrection(stimOnTL, bTLtoMaster);
    
    % call x-file to create stimuli
    SS = stimFile(myScreenInfo, pars.pars);
    stimFrames = cat(3, SS.ImageTextures{:});
    
    framesPerImage = pars.pars(6,1);
    frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / myScreenInfo.FrameRate;
    ft = (frameTimes + noiseOn)';
    
    rows = size(stimFrames,1);
    cols = size(stimFrames,2);
    
    % white squares
    ind = find(stimFrames == 1);
    t = ceil(ind / (rows * cols));
    ind = mod(ind, (rows * cols));
    ind(ind == 0) = rows * cols;
    x_wh = ceil(ind / rows);
    y_wh = mod(ind, rows);
    y_wh(y_wh == 0) = rows;
    time_wh = ft(t,:);
    % black squares
    ind = find(stimFrames == -1);
    t = ceil(ind / (rows * cols));
    ind = mod(ind, (rows * cols));
    ind(ind == 0) = rows * cols;
    x_bl = ceil(ind / rows);
    y_bl = mod(ind, rows);
    y_bl(y_bl == 0) = rows;
    time_bl = ft(t,:);
    
    stimPositions = cell(1,3);
    stimTimes = cell(1,3);
    stimPositions{2} = [repmat(y_wh, 3, 1) repmat(x_wh, 3, 1)];
    stimTimes{2} = reshape(time_wh, [], 1);
    
    stimPositions{3} = [repmat(y_bl, 3, 1) repmat(x_bl, 3, 1)];
    stimTimes{3} = reshape(time_bl, [], 1);
    
    stimTimes{1} = [stimTimes{2}; stimTimes{3}];
    [stimTimes{1}, order] = sort(stimTimes{1}, 'ascend');
    stimPositions{1} = [stimPositions{2}; stimPositions{3}];
    stimPositions{1} = stimPositions{1}(order,:);
    noisePosition = pars.pars(2:5);
    
    stimPseudoTimes = cell(nPseudo,3);
    shifts = randperm(length(frameTimes)-1, nPseudo);
    for j = 1:nPseudo
        ind = find(stimFrames == 1);
        t = ceil(ind / (rows * cols));
        t = mod(t+shifts(j)-1, length(frameTimes))+1;
        stimPseudoTimes{j,2} = reshape(ft(t,:),[],1);
        ind = find(stimFrames == -1);
        t = ceil(ind / (rows * cols));
        t = mod(t+shifts(j)-1, length(frameTimes))+1;
        stimPseudoTimes{j,3} = reshape(ft(t,:),[],1);
        stimPseudoTimes{j,1} = [stimPseudoTimes{j,2}; stimPseudoTimes{j,3}];
    end
    
    %% Load and prepare data for tuning curves
    binSizeGrating = 0.01;
    
    data = load(fullfile(protocolFolder, db(k).subject, db(k).date, ...
        num2str(db(k).expOri), sprintf('%s_%d_%s_parameters.mat', ...
        db(k).date, db(k).expOri, db(k).subject)));
    pars = data.parameters.Protocol;
    
    bTLtoMaster = readNPY(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, sp.name)));
    stimOnTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_onsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
    stimOffTL = readNPY(fullfile(alignDir, ...
        sprintf('mpep_%d_offsets_in_timeline_%d.npy', db(k).expOri, TLexp)));
    stimOn = applyCorrection(stimOnTL, bTLtoMaster);
    stimOff = applyCorrection(stimOffTL, bTLtoMaster);
    stimDur = mean(stimOff - stimOn);
    window = [-0.2 stimDur+0.2];
    binBorders = window(1) : binSizeGrating : window(2);
    binsGrating = binBorders(1:end-1) + binSizeGrating/2;
    stimIDs = repmat((1:pars.npfilestimuli)',pars.nrepeats,1);
    [~,order] = sort(pars.seqnums(:));
    stimSeq = stimIDs(order);
    stimBins = binsGrating>0 & binsGrating<stimDur;
    blank = pars.pars(15,:) == 0;
    directions = pars.pars(6,:);
    directions(blank) = [];
    
    %% Load and prepare data for running correlation
    rotary = double(tl.rawDAQData(:, strcmp({tl.hw.inputs.name}, ...
        'rotaryEncoder')));
    tlTime = applyCorrection(tl.rawDAQTimestamps, bTLtoMaster);
    runSpeed = nonVis.getRunningSpeed_wheel(rotary, tlTime, runningSigma);
    runningTime = runSpeed.t;
    runSpeed = runSpeed.total;
    cmPerUnit = 2*pi * 8.75 / (4 * 1024);
    runSpeed = runSpeed * cmPerUnit;
    binSizeRun = 0.001;
    numD = size(db(k).darkTime,1);
    binsRun = cell(1, numD);
    time = cell(1, numD);
    running = cell(1, numD);
    for d = 1:numD
        binsRun{d} = db(k).darkTime(d,1) : binSizeRun : db(k).darkTime(d,2);
        time{d} = binsRun{d}(1:end-1) + binSizeRun/2;
        running{d} = interp1(runningTime, runSpeed, time{d}, 'pchip');
        time{d}(end+1) = NaN;
        if d>1
            time{d} = time{d} - (time{d}(1)-time{d-1}(end-1)) + 100;
        end
    end
    time = cat(2, time{:});
    time(end) = [];
    
    %% Make plot for whole dataset (recording depth vs. firing rate/visual drive)
    
    folder = fullfile(plotFolder, 'plots_basicFeatures_RFcorrected', ...
        [db(k).subject '_' db(k).date]);
    if ~isdir(folder)
        mkdir(folder)
    end
    
    minAmplitude = 40;
    minSpikes = 100;
    % (1+2) Plot power of spiking in stimulus refresh frequency + mean
    % firing rate
    units = unique(sp.clu);
    stimPowers = NaN(length(units), 1);
    meanRate = NaN(length(units), 1);
    amplitude = NaN(length(units), 1);
    depths = NaN(length(units), 1);
    totalT = sp.st(end)-sp.st(1);
    fs = 1000;
    t = noiseOn + (frameTimes(1) : 1/fs : frameTimes(end));
    if mod(size(t,2), 2) > 0
        t(:,end+1) = t(:,end)+1/fs;
    end
    fr = (0:(size(t,2)/2))./size(t,2).*fs;
    frequencies = myScreenInfo.FrameRate / framesPerImage + [0 -.5 .5];
    freqInds = zeros(size(frequencies));
    for f = 1:length(freqInds)
        freqInds(f) = find(fr >= frequencies(f), 1);
        if diff(abs(fr(freqInds(f)-[1 0]) - frequencies(f))) > 0
            freqInds(f) = freqInds(f)-1;
        end
    end
    for n = 1:length(units)
        ind = sp.clu == units(n);
        amplitude(n) = mean(sp.spikeAmps(ind));
        if amplitude(n) >= minAmplitude || ismember(units(n), goodUnits)
            meanRate(n) = sum(ind) / totalT;
        end
        
        ind = ind & sp.st > t(1) & sp.st < t(end);
        if (sum(ind) < minSpikes || amplitude(n) < minAmplitude) ...
                && ~ismember(units(n), goodUnits)
            continue
        end
        pow = NaN(length(noiseOn), size(t,2)/2+1);
        for trial = 1:length(noiseOn)
            ind = sp.clu == units(n) & sp.st > t(trial,1) & sp.st < t(trial,end);
            s = hist(sp.st(ind), t(trial,:));
            Y = fft(s - mean(s));
            P = abs(Y/size(t,2));
            P = P(1:size(t,2)/2+1);
            P(2:end-1) = 2*P(2:end-1);
            pow(trial,:) = P;
        end
        pow = mean(pow,1);
        stimPowers(n) = pow(freqInds(1)) / mean(pow(freqInds(2:end)));
        depths(n) = mean(sp.spikeDepths(sp.clu == units(n)));
    end
    
    gu = ismember(units, goodUnits);
    
    cols = lines(3);
    ax = zeros(1,3);
    markSz = 15;
    
    figure('Position', [1150 42 720 1074])
    subplot(1,3,1)
    plot(meanRate, depths, 'k.')
    hold on
    plot(meanRate(gu)', depths(gu)', '.', 'Color', cols(1,:), 'MarkerSize', markSz)
    title('Firing rate')
    xlabel('spikes/s')
    ylabel('Depth (um)')
    set(gca, 'XGrid', 'on', 'YGrid', 'on')
    ax(1) = gca;
    
    subplot(1,3,2)
    plot(stimPowers, depths, 'k.')
    hold on
    plot(stimPowers(gu), depths(gu), '.', 'Color', cols(2,:), 'MarkerSize', markSz)
    title('Visual drive')
    xlabel('stim. / neighb. freq.s')
    set(gca, 'XGrid', 'on', 'YGrid', 'on')
    ax(2) = gca;
    
    subplot(1,3,3)
    plot(amplitude, depths, 'k.')
    hold on
    plot(amplitude(gu), depths(gu), '.', 'Color', cols(3,:), 'MarkerSize', markSz)
    title('Spike ampl.')
    xlabel('uV')
    set(gca, 'XGrid', 'on', 'YGrid', 'on')
    ax(3) = gca;
    
    linkaxes(ax, 'y')
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(fullfile(folder, 'depthPlot.jpg'), '-djpeg','-r0')
    
    ylim([min(depths(gu))-20 max(depths(gu))+20])
    
    print(fullfile(folder, 'depthPlot_zoomIn.jpg'), '-djpeg','-r0')
    
    
    close(fig)
    
    %% Make plot for each neuron
    plotRows = 4;
    plotCols = 9;
    plotsTotal = plotRows * plotCols;
    
    for iCell = 1:length(goodUnits)
        % get spiketimes for this cell
        st = sp.st(sp.clu == goodUnits(iCell));
        
        figure('Position', [1 41 1920 1083])
        
        % plot waveform of this cell
        subplot(plotRows, plotCols, 1:plotCols:plotsTotal)
        set(gca, 'Position', [.065 .06 .07 .9])
        hold on
        thisTemp = squeeze(sp.tempsUnW(templates(goodUnits(iCell)+1)+1,:,:));
        [~,peakChan] = max(max(abs(thisTemp),[],1),[],2);
        plotChans = max(1,peakChan-15) : min(length(sp.xcoords),peakChan+15);
        arrayfun(@(x)plot(sp.xcoords(x)+0.3*(1:size(thisTemp,1))', ...
            sp.ycoords(x)+5*thisTemp(:,x), 'k'), plotChans)
        ylabel('Depth (um)')
        title('Waveform')
        
        % plot orientation tuning
        psth = NaN(pars.npfilestimuli, length(binsGrating));
        for s = 1:pars.npfilestimuli
            ind = find(stimSeq == s);
            psth(s,:) = psthAndBA(st, stimOn(ind), window, binSizeGrating);
        end
        meanResp = mean(psth(:,stimBins),2);
        subplot(plotRows,plotCols,reshape([plotCols;2*plotCols]-[1 0],1,[]))
        set(gca, 'Position', [.776 .715 .156 .245])
        imagesc(binsGrating([1 end]),[1 pars.npfilestimuli],psth)
        hold on
        plot([0 0], [0.5 pars.npfilestimuli+.5], 'w:', 'LineWidth', 2)
        plot([1 1].*stimDur, [0.5 pars.npfilestimuli+.5], 'w:', 'LineWidth', 2)
        xlabel('Time from stim onset (s)')
        ylabel('Stimulus')
        title('Direction tuning')
        colormap hot
        colorbar('Position', [.941 .715 .014 .245])
        subplot(plotRows,plotCols,reshape([3;4].*plotCols-[1 0],1,[]))
        set(gca, 'Position', [.776 .395 .156 .245])
        plot(directions, meanResp(~blank), 'o-k', 'MarkerFaceColor', 'k')
        hold on
        plot(directions([1 end]), [1 1] .* meanResp(blank), 'k:')
        xlim([directions(1)-10 directions(end)+10])
        xlabel('Direction')
        ylabel('Firing rate')
        set(gca,'box','off')
        
        % plot running correlation
%         spikerate = cell(1, numD);
%         rho = NaN(1, numD);
%         pVal = NaN(1, numD);
%         sigma = round(spikingSigma / binSizeRun);
%         run = running;
%         for d = 1:numD
%             stRun = st(st >= binsRun{d}(1) & st <= binsRun{d}(end));
%             spikerate{d} = histcounts(stRun, binsRun{d}) ./ binSizeRun;
%             win = normpdf(-5*sigma : 5*sigma, 0, sigma);
%             spikerate{d} = conv(spikerate{d}, win, 'same');
%             [rho(d),pVal(d)] = corr(run{d}', spikerate{d}');
%             run{d}(end+1) = NaN;
%             spikerate{d}(end+1) = NaN;
%         end
%         run = cat(2, run{:});
%         run(end) = [];
%         spikerate = cat(2, spikerate{:});
%         spikerate(end) = [];
%         
%         ax = zeros(1,2);
%         subplot(plotRows,plotCols,plotsTotal+(-plotCols+2:0))
%         set(gca, 'Position', [.218 .06 .714 .103])
%         plot(time, spikerate, 'k')
%         ax(2) = gca;
%         ylabel('Spikes/s')
%         set(gca,'box','off')
%         xlabel('Time (s)')
%         subplot(plotRows,plotCols,plotsTotal-plotCols+(-plotCols+2:0))
%         set(gca, 'Position', [.218 .204 .714 .103])
%         plot(time, run)
%         ax(1) = gca;
%         set(gca,'box','off')
%         linkaxes(ax, 'x')
%         xlim(time([1 end]))
%         ylabel('Running speed (cm/s)')
%         title(['Correlation with running (r=' sprintf('%.2f ',rho) ', p=' ...
%             sprintf('%.3f ',pVal) ')'])
        
        % plot RF of this cell
        spatialRFs = cell(length(RFtypes),1);
        timeCourses = cell(length(RFtypes),1);
        residuals = cell(length(RFtypes),1);
%         [RF, RFtime] = whiteNoise.getRF_spikes(st, stimFrames, frameTimes, ...
%             noiseOn, [0 0], RFtypes);
%         spatialRFs = RF';

        psths = cell(length(RFtypes),1);
        RF = cell(length(RFtypes),1);
        for type = 1:length(RFtypes)
            [~, stats, psths{type}] = sparseNoiseRF(st, stimTimes{type}, ...
                stimPositions{type}, params);
            
            pseudoPsths = cell(1,nPseudo);
            for j = 1:nPseudo
                [~,~, pseudoPsths{j}] = sparseNoiseRF(st, stimPseudoTimes{j,type}, ...
                    stimPositions{type}, params);
            end
            
            RF{type} = psths{type} - mean(cat(3, pseudoPsths{:}), 3);
        end
        RFtime = stats.timeBins;
        binSize = mean(diff(RFtime));
        
        for type = 1:length(RFtypes)
            [U,S,V] = svd(RF{type} - nanmean(RF{type}(:)), 'econ');
            diffs = -diff(diag(S));
            [d,numComp] = max(diffs(1:2));
            if d < 4*median(diffs)
                numComp = 1;
            end
            for comp = 1:numComp
                spatialRFs{type,comp} = reshape(U(:,comp), size(stimFrames,1), ...
                    size(stimFrames,2)) .* S(comp,comp);
                timeCourses{type,comp} = V(:,comp);
                [~,peakInd] = max(abs(U(:,comp)));
                if U(peakInd,comp) < 0
                    spatialRFs{type,comp} = -spatialRFs{type,comp};
                    timeCourses{type,comp} = -timeCourses{type,comp};
                end
                if comp == 1
                    residuals{type,comp} = RF{type} - (U(:,comp) .* ...
                        S(comp,comp) * V(:,comp)' + nanmean(RF{type}(:)));
                else
                    residuals{type,comp} = residuals{type,comp-1} - ...
                        (U(:,comp) .* S(comp,comp) * V(:,comp)');
                end
            end
        end
        
        legs = {'comp 1', 'comp 2'};
        subplot(plotRows,plotCols,2:3)
        set(gca, 'Position', [.218 .86 .156 .1])
        hold on
        for comp = 1:size(spatialRFs,2)
            if ~isempty(timeCourses{1,comp})
                plot(RFtime, timeCourses{1,comp}, 'Color', ...
                    colors(type,:))
            end
        end
        axis tight
        set(gca, 'box', 'off')
        xlabel('Time from stim onset (s)')
        title('RF time course')
        legend(legs(1:size(spatialRFs,2)), ...
            'Position', [.38 .86 .046 .1]);
        
        m = max(abs([RF{1}(:); residuals{1,1}(:)]));
        subplot(plotRows, plotCols, plotCols+(4:7))
        imagesc([1 size(RF{1},1)], [RFtime(1)-binSize/2 ...
            RFtime(end)+binSize/2], RF{1}', [-m m])
        colorbar
        title('PSTH')
        ylabel('Time')
        colormap(gca,cm)
        
        for comp = 1:size(spatialRFs,2)
            if isempty(spatialRFs{1,comp})
                continue
            end
            subplot(plotRows,plotCols,plotCols*(1+comp)+[2 3])
            m = max(abs(spatialRFs{1,comp}(:)));
            imagesc(noisePosition([1 2]), noisePosition([3 4]), ...
                spatialRFs{1,comp}, [-m m])
            colorbar
            if comp==1
                title('Spatial RF')
            end
            ylabel(['Comp ' num2str(comp)], 'FontSize', 13, 'FontWeight', 'bold')
            colormap(gca,cm)
            subplot(plotRows,plotCols,plotCols*(1+comp)+[4 7])
            m = max(abs([RF{1}(:); residuals{1,1}(:)]));
            imagesc([1 size(RF{1},1)], [RFtime(1)-binSize/2 ...
                RFtime(end)+binSize/2], residuals{1,comp}', [-m m])
            colorbar
            if comp==1
                title('Residuals')
            end
            xlabel('Pixel')
            colormap(gca,cm)
        end
            
        annotation('textbox', [.47 .95 .17 .04], 'String', sprintf(...
            'Neuron %d (good unit %d)', goodUnits(iCell), iCell), ...
            'FontSize', 15, 'FontWeight', 'bold', ...
            'LineStyle', 'none', 'HorizontalAlignment', 'left')
        ind = find(units == goodUnits(iCell));
        annotation('textbox', [.47 .88 .17 .07], 'String', sprintf(...
            'Depth: %d um\nVisual drive: %.2fx\nSpike ampl.: %d uV', ...
            round(depths(ind)), stimPowers(ind), round(amplitude(ind))), ...
            'FontSize', 11, ...
            'LineStyle', 'none', 'HorizontalAlignment', 'left')
        
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(fullfile(folder, sprintf('neuron%03d.jpg', ...
            goodUnits(iCell))), '-djpeg','-r0')
        close(fig)
    end        
    
%     goodSpikes = ismember(sp.clu, goodUnits);
%     psthViewer(sp.st(goodSpikes), sp.clu(goodSpikes), stimOn, window, ...
%         stimSeq)
end