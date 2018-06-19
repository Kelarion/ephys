

% New RF online, using signals RF mapper and LFP rather than mpep and
% spikes... 

%% paths
addpath(genpath('C:\Users\Experiment\Documents\GitHub\spikes'))
addpath(genpath('C:\Users\Experiment\Documents\GitHub\Rigbox'))


%% specify files, other parameters

% 
% localDataDir = 'F:\data';
% mouseName = 'SS093';
% thisDate = datestr(now, 'yyyy-mm-dd');
% tag = 'K1';
% expNum = 1;
% localDataDir = 'G:\data';
% mouseName = 'SS093';
% thisDate = datestr(now, 'yyyy-mm-dd');
% tag = 'K2';
% expNum = 1;
localDataDir = 'J:\data';
mouseName = 'SS093';
thisDate = datestr(now, 'yyyy-mm-dd');
tag = 'K3';
expNum = 1;

FsOrig = 2500;
downSampFactor = 5;
Fs = FsOrig/downSampFactor;

exclChans = [37    76   113   152   189   228   265   304   341   380]; %internal refs
q = 1:384;
inclChans = q(~ismember(q, exclChans));
useChans = inclChans(1:8:end);

% inclSamps = 2.85e6:4.4e6; % this is just for my test file


computeWin = [-0.05 0.2]; % window around stimulus events to compute average LFP


%% load LFP data
fprintf(1, 'loading data\n');
tic

% lfD = dir(fullfile(localDataDir, mouseName, thisDate, tag, '*.lf.bin'));
% lfFn = fullfile(localDataDir, mouseName, thisDate, tag, lfD.name);
lfD = dir(fullfile(localDataDir, mouseName, thisDate, '*.lf.bin'));
lfFn = fullfile(localDataDir, mouseName, thisDate, lfD.name);

nChans = 385;
nSampToRead = floor(lfD.bytes/2/nChans);

fprintf(1, 'file appears to have %.2f sec (%.2f min) of data in it\n', nSampToRead/Fs, nSampToRead/Fs/60);

fid = fopen(lfFn);
% could improve this by skipping the channels I will drop anyway, perhaps
% by memory mapping. would reduce memory usage, probably not increase speed
rawData = fread(fid, [nChans nSampToRead], '*int16');
fclose(fid);

syncDat = rawData(end,:);
lfDat = rawData(useChans,:);

% syncDat = syncDat(inclSamps);
% lfDat = lfDat(:,inclSamps);

syncDat = syncDat(1:downSampFactor:end);
lfDat = lfDat(:,1:downSampFactor:end);

clear rawData

toc; tic
fprintf(1, 'preprocess data\n');

% subtract channel medians because data comes offset. 
lfDat = bsxfun(@minus, lfDat, median(lfDat,2));

%% extract times of syncEvents
% ds = diff(syncDat);
% syncUp = find(ds==1);
% syncDown = find(ds==-1);
% syncEvents = sort([syncUp syncDown])/Fs;
% syncEvents = syncEvents(:); % make column
% syncEvents = [0; syncEvents]; 
eventTimes = spikeGLXdigitalParse(syncDat, Fs);
syncEvents = [0;eventTimes{1}{1}];

t = (0:size(lfDat,2)-1)/Fs;

toc
%% load stimulus information
fprintf(1, 'load stimulus info\n');tic
load(dat.expFilePath(mouseName, thisDate, expNum, 'block', 'master'))

sw = block.stimWindowUpdateTimes;
sw = sw(:); % make column

[stimTimeInds, stimPositions, stimArray] = ...
    computeSparseNoiseSignals(block);

% these come with two cells for the pixel turning on or turning off. Here I
% assume they only turn on (black and white version of noise).
stimTimeInds = stimTimeInds{1}; 
stimPositions = stimPositions{1};
toc
%% work out alignment between block and FPGA
% Method is to take the same number of FPGA events as we need to match the
% number from the block, and work along until we find ones that fit. "Fit"
% means the difference between the diff's of the two cannot exceed one
% frame, i.e. they move step-for-step together. No frames can be missed, or
% this would fail. 
fprintf(1, 'find alignment\n');tic
assert(numel(syncEvents)>=numel(sw), ...
    'not enough sync events in the fpga. something wrong! cannot proceed.\n');
    
startInd = 1;
match = false;
while ~match && startInd+numel(sw)-1<=numel(syncEvents)
    theseSE = syncEvents(startInd:startInd+numel(sw)-1);
    
    match = sum((diff(sw)-diff(theseSE))>0.02)<10; % the diffs have to be within one frame (20ms) but allow a few to be off for random reasons
    startInd = startInd+1;
end

assert(match, 'could not find matching events in sync and block\n');

stimTimes = theseSE(stimTimeInds);
toc
%% finish computing stimulus information given alignment

[stimTimes, ii] = sort(stimTimes);
stimPositions = stimPositions(ii,:);

xPos = unique(stimPositions(:,1)); nX = length(xPos);
yPos = unique(stimPositions(:,2)); nY = length(yPos);

%% compute the RF for each LFP trace given times

fprintf(1, 'computing stim-locked responses\n'); tic

nCh = size(lfDat,1); 
winSamps = computeWin(1):1/Fs:computeWin(2);
nWS = length(winSamps);
allResp = zeros(nCh, nX, nY, nWS);

periEventTimes = bsxfun(@plus, stimTimes, winSamps);
theseLockedAll = interp1(t, double(lfDat)', periEventTimes);

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

toc

%% svds

% fprintf(1, 'computing maps using svd\n'); tic
% rfMaps = zeros(nCh, nX, nY);
% timeCourses = zeros(nCh, nWS);
% 
% for c = 1:nCh
%     thisResp = reshape(allResp(c,:,:,:), nX*nY, nWS);
%     bsl = mean(thisResp(:,1)); % take the first bin as the baseline
%     [U,S,V] = svd(thisResp - bsl,'econ');
%     rfMapVec = U(:,1);
%     rfMaps(c,:,:) = reshape(rfMapVec,nX, nY);
%     timeCourses(c,:) = V(:,1)';
% %     Scalar = S(1,1);
% %     Model = rfMapVec*timeCourse*Scalar + bsl;
% %     Residual = thisResp - Model;
% end
% 
% toc

%% using peaks

fprintf(1, 'computing maps using peaks\n'); tic
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

toc

%% plot responses
fprintf(1, 'plotting maps\n');tic
f = figure;
for c = 1:nCh
    cPlotInd = nCh-c+1;
    subtightplot(nCh, 2, (cPlotInd-1)*2+1, 0.001, 0.001, 0.001);
    imagesc(yPos, xPos, squeeze(rfMaps(c,:,:)));
    hold on; plot(90,0, 'ro');
    hold on; plot(-90,0, 'ro');
    axis image; axis off
    cax = caxis();
    caxis([-max(abs(cax)) max(abs(cax))]);
    %ylabel('c')
    
%     subtightplot(nCh, 2, (cPlotInd-1)*2+2, 0.01, 0.01, 0.01);
%     plot(winSamps, timeCourses(c,:));
%     ylim([min(timeCourses(c,:)) max(timeCourses(c,:))]);
%     hold on;
%     plot([0 0], ylim(), 'k--');
%     xlim(computeWin);
%     axis off;
%     drawnow;
    
end

subtightplot(1,2,2, 0.001, 0.001, 0.001);
imagesc(winSamps, 1:nCh, timeCourses);
xlabel('time (s)')
hold on; plot([0 0], [0 nCh+1], 'k--');
axis off
set(gca, 'YDir', 'normal')
cax = caxis();
caxis([-max(abs(cax)) max(abs(cax))]);
% colormap(colormap_RedWhiteBlue)
toc

%% plot a single channel

c = 33; % likely in cortex

figure; 
subplot(2,1,1);
imagesc(yPos, xPos, squeeze(rfMaps(c,:,:)));
hold on; plot(90,0, 'go');
axis image; %axis off
cax = caxis();
caxis([-max(abs(cax)) max(abs(cax))]);
% colormap(colormap_RedWhiteBlue)

subplot(2,1,2);
plot(winSamps, timeCourses(c,:));
xlabel('time (s)')
hold on; plot([0 0], ylim(), 'k--');
xlim(computeWin)
box off