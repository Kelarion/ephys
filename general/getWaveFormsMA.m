function wf = getWaveFormsMA(gwfparams)
% function wf = getWaveForms(gwfparams)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% Contributed by C. Schoonover and A. Fink
% Chopped, not slopped, by M. Alleman
%
% % EXAMPLE INPUT
% gwfparams.dataDir = '/path/to/data/';    % Folder with actual data files
% gwfparams.ksDir = '/path/to/chanmap/';   % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 385;                     % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-31 35];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 1000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
% gwfparams.nMaxCh = [7 8];                % Number of channels to save, before and after max-amplitude channel
% gwfparams.maxSite = [3,7,19, ... ];      % Max-amplitude sites
%
% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nMaxCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nMaxCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);             
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
% nSamp = double(max(gwfparams.spikeTimes) + gwfparams.wfWin(2));
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'}); 
chMap = readNPY([gwfparams.ksDir 'channel_map.npy'])+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChSave = sum(gwfparams.nMaxCh)-1;

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
numUnits = size(unitIDs,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,nChSave,wfNSamples);
waveFormsMean = nan(numUnits,nChSave,wfNSamples);
tic;
for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    maxChan = gwfparams.maxSite(curUnitInd); % max-amplitude channels
%     if maxChan>gwfparams.nMaxCh(1) && maxChan<(length(chMap)-gwfparams.nMaxCh(2))
%         mySites = maxChan-gwfparams.nMaxCh(1):maxChan+gwfparams.nMaxCh(2);
%     elseif maxChan <= gwfparams.nMaxCh(1) % if max site is at the beginning
%         mySites = 1:maxChan+gwfparams.nMaxCh(2)+;
%     else % if it's at the end of the probe
%         mySites = maxChan-gwfparams.nMaxCh(1):length(chMap);
%     end
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimesRP(spikeTimesRP > nSamp - gwfparams.wfWin(end)) = [];
    spikeTimesRP(spikeTimesRP < gwfparams.wfWin(1)) = [];
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
        tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap(maxChan),:);
    end
    waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:)));
%     disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) ' in ' num2str(toc) ' seconds.']);
end

% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;

end
