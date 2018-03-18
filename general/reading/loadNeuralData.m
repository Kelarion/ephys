function spks = loadNeuralData(ksDir,dataDir)
% spks = loadNeuralData(ksDir[,dataDir])
%
% I think this loads basically the same information as loadAllKsDir, but on
% a single-tag basis and somehow is a bit more flexible for my purposes.
%
% if dataFolder isn't provided, assumes is the parent directory of ksDir

if nargin < 2
    parts = split(ksDir,filesep);
    dataDir = join(parts(1:end-2),filesep);
    dataDir = dataDir{:};
end

chanMap = readNPY([ksDir 'channel_map.npy']);
spks = loadKSdir(ksDir); % loads all sorts of things (check out the function)
spks.channel_map = double(chanMap);

filenamestruct = dir([dataDir '\' spks.dat_path]);
if isempty(filenamestruct)
    error('No data in provided directory')
end
dataTypeNBytes = numel(typecast(cast(0, spks.dtype), 'uint8')); % determine number of bytes per sample
nSampDat = filenamestruct.bytes/(spks.n_channels_dat*dataTypeNBytes); % to get number of samples

foo = split(spks.dat_path,'.ap');
lfpName = [foo{1} '.lf.bin'];
lfpnamestruct = dir([dataDir '\' lfpName]);
dataTypeNBytes = numel(typecast(cast(0, spks.dtype), 'uint8')); % determine number of bytes per sample
nSampLFP = lfpnamestruct.bytes/(spks.n_channels_dat*dataTypeNBytes);
lfpFs = round(spks.sample_rate/(nSampDat/nSampLFP)); % lfp is downsampled by a certain factor
                                                 % so this should already be integer anyway
spks.nSampDat = nSampDat;
spks.lfp_path = lfpName;
spks.nSampLFP = nSampLFP;
spks.lfp_sample_rate = lfpFs; 

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
