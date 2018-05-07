function spks = loadNeuralData(ksDir,dataDir)
% spks = loadNeuralData(ksDir[,dataDir])
%
% I think this loads basically the same information as loadAllKsDir, but on
% a single-tag basis and with some other info that I like for some reason.
%
% if dataDir isn't provided, assumes is the parent directory of ksDir
%
% Specifically, what this does in addition to loadAllKsDir is:
%
% nSampDat:         number of 
% lfp_path:         just the path to the lfp file
% nSampLFP:         same deal as above
% lfp_sample_rate:  probably not super useful
% cluTemps:         output of findTempForEachClu [nClu×1] 
% isNeuron:         spikes coming from non-noise clusters [nSpk×1 logical]
% isGood:           spike coming from good clusters [nSpk×1 logical]

if nargin < 2
    parts = split(ksDir,filesep);
    dataDir = join(parts(1:end-2),filesep);
    dataDir = dataDir{:};
end

chanMap = readNPY([ksDir 'channel_map.npy']);
spks = loadKSdir(ksDir); % loads all sorts of things (check out the function)
if ~isfield(spks,'gain')
    spks.gain = 2.3438;
    spks.gainLFP = 4.6875;
end
spks.channel_map = double(chanMap);

filenamestruct = dir([dataDir '\' spks.dat_path]);
if isempty(filenamestruct)
    error('No data in provided directory')
end
dataTypeNBytes = numel(typecast(cast(0, spks.dtype), 'uint8')); % determine number of bytes per sample
nSampDat = filenamestruct.bytes/(spks.n_channels_dat*dataTypeNBytes); % to get number of samples

foo = split(spks.dat_path,'.ap');
lfpName = [foo{1} '.lf.bin'];
if exist([dataDir lfpName],'file')
    lfpnamestruct = dir([dataDir '\' lfpName]);
    dataTypeNBytes = numel(typecast(cast(0, spks.dtype), 'uint8')); % determine number of bytes per sample
    nSampLFP = lfpnamestruct.bytes/(spks.n_channels_dat*dataTypeNBytes);
    lfpFs = round(spks.sample_rate/(nSampDat/nSampLFP)); % lfp is downsampled by a certain factor
    % so this should already be integer anyway
    spks.lfp_path = lfpName;
    spks.nSampLFP = nSampLFP;
    spks.lfp_sample_rate = lfpFs;
end
spks.nSampDat = nSampDat; 

spks.cluTemps = findTempForEachClu(double(spks.clu),double(spks.spikeTemplates));

[spks.spikeAmps, spks.spikeDepths, spks.templateDepths, ~, spks.tempsUnW, ~, ~] = ...
    templatePositionsAmplitudes(spks.temps, spks.winv, spks.ycoords,  ...
    spks.spikeTemplates, spks.tempScalingAmps);

if any(spks.cgs==2)
    spks.isNeuron = ismember(spks.clu,spks.cids(spks.cgs==2 | spks.cgs==1));
    spks.isGood = ismember(spks.clu,spks.cids(spks.cgs==2));
else
    spks.isNeuron = ismember(spks.clu,spks.cids);
    spks.isGood = [];
end
