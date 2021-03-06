

function [rfMap, stats, allPSTH] = sparseNoiseRF_MA(spikeTimes, stimTimes, stimPositions, params)
% function [rfmap, stats] = sparseNoiseRF(spikeTimes, stimTimes, stimPositions, params)
%
% Assumes that stimPositions describe a rectangle and are evenly spaced.
%
% - spikeTimes is nSpikes x 1
% - stimTimes is nStimulusEvents x 1
% - stimPositions is nStimulusEvents x 2 (x and y coordinates)
%
% params can have: 
% - makePlots - logical
% - useSVD - logical, whether to compute the RF by using SVD on the PSTHs for all
%   stimulus responses. If not, will just count spikes in a window. 
% - countWindow - 1x2, start and end time relative to stimulus onset to
%   consider
% - binSize - 1x1, to use for making rasters
% - fit2Dgauss - logical, try to fit a 2D gaussian
% - nShuffle - scalar, number of permutations for generating null
%   distribution of RF peaks (default empty --> no perm. test, returns
%   peak z-score relative to RF rather than null).
%
% TODO
% - show some rasters of the peak and of an average site
% - implement statistical test(s)
% - make colormap divergent so suppression is clear
% - show also raw count map even when using model

% default parameters
p.makePlots = true; 
p.useSVD = true;
p.countWindow = [0 0.2];
p.binSize = 0.01;
p.fit2Dgauss = false;
p.nShuffle = [];

stats = [];

fn = fieldnames(p);
for f = 1:length(fn)
    if isfield(params, fn{f}) && ~isempty(params.(fn{f}))
        p.(fn{f}) = params.(fn{f});
    end
end

useShuffle = ~isempty(p.nShuffle);

[stimTimes, ii] = sort(stimTimes);
stimPositions = stimPositions(ii,:);

xPos = unique(stimPositions(:,1)); nX = length(xPos);
yPos = unique(stimPositions(:,2)); nY = length(yPos);

if p.useSVD
    timeBins = [p.countWindow(1):p.binSize:p.countWindow(2)];
    timeBins = timeBins(1:end-1)+p.binSize/2;
    nBins = length(timeBins);
    thisRF = zeros(nX, nY, nBins);
    stats.timeBins = timeBins;
else
    timeBins = [p.countWindow(1):p.binSize:p.countWindow(2)];
    timeBins = timeBins(1:end-1)+p.binSize/2;
    nBins = length(timeBins);
    thisRF = zeros(nX,nY);
    psth = zeros(nX, nY, nBins);
    stats.timeBins = timeBins;
end

[~, ~, ~, ~, spikeCounts, ba] = psthAndBA(spikeTimes, stimTimes, p.countWindow, p.binSize);

for x = 1:nX
    for y = 1:nY
        theseStims = stimPositions(:,1)==xPos(x) & stimPositions(:,2)==yPos(y);
        if sum(theseStims)==0
            warning('sparseNoiseRF: stimulus at %d, %d was never shown', xPos(x), yPos(y));
        else
            if p.useSVD
                thisRF(x,y,:) = mean(ba(theseStims,:));
            else
                thisRF(x,y) = mean(spikeCounts(theseStims));
                psth(x,y,:) = mean(ba(theseStims,:));
            end
        end
    end
end





if p.useSVD
    allPSTH = reshape(thisRF, nX*nY, size(thisRF,3));
    bsl = mean(allPSTH(:,1)); % take the first bin as the baseline
    [U,S,V] = svd(allPSTH - bsl,'econ');
    
    rfMapVec = U(:,1);
    timeCourse = V(:,1)';
    
    % here I try to flip the SVD the correct way. For a neuron that's
    % activated, this is easy: we want the peak space/time point to be two
    % positive numbers rather than two negative ones. For a neuron that's
    % suppressed it's more complicated because it's reasonable to describe
    % that neuron either with a positive time course and negative field, or
    % negative time course and positive field. Here I make the decision
    % to go with the former, since I mostly just look at the maps later. So
    % if the peak point in the model is negative in time but positive in 
    % space, we will flip it.     
    %     +time, neg peak = +/-, no flip
    %     +time, pos peak = +/+, no flip
    %     -time, neg peak = -/+, flip
    %     -time, pos peak = -/-, flip
    % So the result of all this is just that we never want the time course
    % to have a negative peak.
    % * What if it has both, at different points in time? 
    peakTimePoint = find(abs(timeCourse)==max(abs(timeCourse)),1);    
    ts = sign(timeCourse(peakTimePoint));
    if ts<0
        rfMapVec = -rfMapVec;
        timeCourse = -timeCourse;
    end
        
    rfMap = reshape(rfMapVec,nX, nY);
    Scalar = S(1,1);
    Model = rfMapVec*timeCourse*Scalar + bsl;
    Residual = allPSTH - Model;
    stats.timeCourse = timeCourse;
    stats.scalar = Scalar;
    stats.bsl = bsl;
else
    rfMap = thisRF;
    allPSTH = reshape(psth, nX*nY, size(psth,3));
    stats.scalar = 1;
    stats.bsl = mean(allPSTH(:,1));
%     stats.timeCourse;
end


% statistical test(s)
% - shuffle test: idea is that if you relabel each stimulus event with a different
% position, on what proportion of relabelings do you get a peak as big as
% the one actually observed? or can calculate this analytically from the
% distribution of all spike counts.
% - simplest: just the z-score of the peak relative to the whole population

if useShuffle
    [~,~,pos_grps] = unique(stimPositions,'rows');
    all_pos = 1:(nX*nY);
    spk_cnts = clusterAverage(pos_grps,spikeCounts);
    maxDev = (max(spk_cnts)-mean(spk_cnts));%./std(spk_cnts);
    minDev = (min(spk_cnts)-mean(spk_cnts));%./std(spk_cnts);
    truePeak =  max(abs([minDev maxDev])); % use the unstandardized max deviation

    nullPeaks = zeros(p.nShuffle,1);
    for iShuf = 1:p.nShuffle
        sigma = randperm(length(stimPositions));
        
        fooCounts = zeros(nX*nY,1);
        for r = 1:(nX*nY)
            theseStims = (pos_grps(sigma) == all_pos(r));
            [n,~] = histdiff(spikeTimes, stimTimes(theseStims),  p.countWindow);
            fooCounts(r) = n./sum(theseStims);
        end
        
        maxDev = (max(fooCounts)-mean(fooCounts));%./std(fooCounts);
        minDev = (min(fooCounts)-mean(fooCounts));%./std(fooCounts);
        nullPeaks(iShuf) =  max(abs([minDev maxDev]));

    end
    stats.peakZscore = truePeak/std(spk_cnts);
    stats.z_peak = (truePeak - mean(nullPeaks))./std(nullPeaks);
else
    maxZ = (max(rfMap(:))-mean(rfMap(:)))./std(rfMap(:));
    minZ = (min(rfMap(:))-mean(rfMap(:)))./std(rfMap(:));
    stats.peakZscore = max(abs([minZ maxZ]));
%     stats.p_zscore = [];
    stats.metaZ = [];
end

if p.fit2Dgauss
    
    if abs(minDev)>maxDev
        % peak is negative - neuron is suppressed        
        useRF = -rfMap; 
    else
        useRF = rfMap;
    end
    useRF = useRF-quantile(useRF(:),0.25);
    
    x = fit2dGaussRF(yPos, xPos, useRF, false);
    stats.fit2D = x;
    
end


if p.makePlots
    
    
    figure; 
    
    
    if p.useSVD
        nCol = 2;
        
        subplot(3,nCol,2);
        plot(timeBins, timeCourse);
        xlim([timeBins(1) timeBins(end)]);
        title('response time course')
        xlabel('time')
        
        ax = subplot(3,nCol,4);
        imagesc(1:size(allPSTH,1), timeBins,  allPSTH');
%         cax = caxis();
%         caxis(max(abs(cax))*[-1 1]);
        colormap(colorcet('L3'));
%         axis image
        title('all PSTHs, ordered');
        xlabel('space')
        ylabel('time')
        
        subplot(3,nCol,6); 
        imagesc(Residual');
%         cax = caxis();
%         caxis(max(abs(cax))*[-1 1]);     
%         axis image
        title('residual')
        xlabel('space')
        ylabel('time')
    else
        nCol = 1;
    end
    
    subplot(3,nCol,1);
    imagesc(yPos, xPos, rfMap);
    axis image
    title(sprintf('map, peakZ = %.2f', stats.peakZscore));
    %colormap hot
%     colormap(colorcet('D9'));
%     cax = caxis();
%     caxis(max(abs(cax))*[-1 1]);
    
    if p.fit2Dgauss
        subplot(3,nCol,(nCol-1)*2+1)
        imagesc(yPos, xPos, rfMap);
        title('map with fit');
        colormap hot
        
        
        %subplot(3,nCol,(nCol-1)*3+1)
        t = linspace(0,2*pi,1000);
        xx = x(2)+x(3)*sin(t+x(6));
        yy = x(4)+x(5)*cos(t);
        
        hold on;
        plot(xx,yy,'g', 'LineWidth', 2.0);
    end
    
    
    
end

