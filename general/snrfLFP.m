function [timeCourses, rfMaps, channels] = snrfLFP(lfpFile,stimPositions,stimTimes,useChans,computeWin)
% [timeCourses, rfMaps, channels] = snrfLFP(lfpFile,stimPosition,stimTimes[,useChans,computeWin])
% 
% lfpFile is a memmap of the LFP file. make sure 'stimTimes' is in the
% timeframe of this specific probe.
%
% modified from rfOnlineSigLF (in from_Nick folder)

if nargin<6 
    computeWin = [-0.05 0.2];
end
nDepthBins = 8;
dsFactor = 1;

lfpFs = 2500/dsFactor;
nChLFP = size(lfpFile.Data.x,1);
xPos = unique(stimPositions(:,1)); nX = length(xPos);
yPos = unique(stimPositions(:,2)); nY = length(yPos);

if nargin<5
    useChans = 1:nChLFP;
else
    nChLFP = length(useChans);
end

channels = useChans(1:nDepthBins:end);
nCh = length(channels);

tsnrf = min(stimTimes)-3:1/lfpFs:max(stimTimes)+3;
whenNoise = round(tsnrf*lfpFs);
winSamps = computeWin(1):1/lfpFs:computeWin(2);
nWS = length(winSamps);
allResp = zeros(nCh, nX, nY, nWS);

periEventTimes = bsxfun(@plus, stimTimes, winSamps);
% B = binmat(nChLFP,nDepthBins)/nDepthBins; lfDat = B*double(lfpFile.Data.x(useChans,whenNoise));
lfDat = double(lfpFile.Data.x(channels,whenNoise));
lfDat = bsxfun(@minus,lfDat,median(lfDat,2)); 
theseLockedAll = interp1(tsnrf, lfDat', periEventTimes);

theseStims = {};
for x = 1:nX
    for y = 1:nY
        theseStims{x,y} = stimPositions(:,1)==xPos(x) & stimPositions(:,2)==yPos(y);
    end
end
    
for c = 1:nCh        
    for x = 1:nX
        for y = 1:nY
            allResp(c,x,y,:) = mean(theseLockedAll(theseStims{x,y},:,c));
        end
    end    
end

rfMaps = zeros(nCh, nX, nY);
timeCourses = zeros(nCh, nWS);

for c = 1:nCh
    thisResp = reshape(allResp(c,:,:,:), nX*nY, nWS);
    bsl = nanmean(thisResp(:,1)); % take the first bin as the baseline
    
    respNorm = thisResp-bsl;
    
    respNorm(isnan(respNorm)) = 0;
    
    % find the peak stimulus and time point
    [iStim,iTimePoint] = find(abs(respNorm)==max(abs(respNorm(:))),1);
    
    % take the map to be the value at that time point for every stimulus
    % position
    rfMaps(c,:,:) = reshape(respNorm(:,iTimePoint),nX, nY);
    
    % take the time course to be the time course of the best stimulus
    timeCourses(c,:) = respNorm(iStim,:);
    
end


end
