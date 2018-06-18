function [stimTimes, stimPositions] = eventTimesLocal(dataDir,blockStruct,Fs)
% eTimes = eventTimesLocal(ephysFile,snBlock)
% 
% Return the sparse noise event times in blockStruct in the time frame of
% the recording specified by ephys_fname.
%
% Modified from rfOnlineSigLF (in from_Nick)

if nargin<3, Fs = 2500; end

ephysD = dir([dataDir '\*lf.bin']);
assert(~isempty(ephysD),'no data file in given directory!')

% ephys_fname = [dataDir ephysD.name];

% nChans = 385;
% nSampToRead = floor(ephysD.bytes/2/nChans);

syncDat = extractSyncChannel(dataDir,385,385);

% ephysFile = memmapfile(ephys_fname, 'Format', {'int16', [nChans nSampToRead], 'x'});
% syncDat = ephysFile.Data.x(end,:);

%% extract times of syncEvents

eventTimes = spikeGLXdigitalParse(syncDat, Fs);
syncEvents = [0;eventTimes{1}{1}];

%% load stimulus information
sw = blockStruct.stimWindowUpdateTimes;
sw = sw(:); % make column

[stimTimeInds, stimPositions] = ...
    computeSparseNoiseSignals(blockStruct);

% these come with two cells for the pixel turning on or turning off. Here I
% assume they only turn on (black and white version of noise).
stimTimeInds = stimTimeInds{1}; 
stimPositions = stimPositions{1};

%% work out alignment between block and FPGA
% Method is to take the same number of FPGA events as we need to match the
% number from the block, and work along until we find ones that fit. "Fit"
% means the difference between the diff's of the two cannot exceed one
% frame, i.e. they move step-for-step together. No frames can be missed, or
% this would fail. 
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

