%% load timeline and blocks
tl = load([infoRoot expFolders num2str(db(k).tlExp) '\' tlInfo '_Timeline.mat']);
blks.tl = tl.Timeline;
cwname = sprintf('%s_%d_%s_Block.mat',db(k).date,db(k).cwExp,db(k).mouse_name);
tmpcw = load([infoRoot expFolders num2str(db(k).cwExp) '\' cwname]);
blks.cw = tmpcw.block;
passivename = sprintf('%s_%d_%s_Block.mat',db(k).date,db(k).passiveExp,db(k).mouse_name);
tmppas = load([infoRoot expFolders num2str(db(k).passiveExp) '\' passivename]);
blks.passive = tmppas.block;

%% sparse noise information
if exist(alfFolder,'dir')
    sn.stimTimes = readNPY([alfFolder 'sparseNoise.times.npy']);
    sn.stimPosition = readNPY([alfFolder 'sparseNoise.positions.npy']);
end

%% load alignment info
alignToRef = dir([alignFolder 'correct_timeline_*.npy']);
[~, refInd] = regexp(alignToRef.name,'ephys_');
aln.refTag = alignToRef.name(refInd+[1:2]);

if strcmp(thisTag,aln.refTag)
    aln.ref2ephys = [1 0];
else
    aln.ref2ephys = readNPY([alignFolder sprintf('correct_ephys_%s_to_ephys_%s.npy', ...
        thisTag,aln.refTag)]); % correct from ref ephys to current ephys
end
aln.tl2ref = readNPY([alignToRef.folder '\' alignToRef.name]);
aln.cw2tl = readNPY([alignFolder sprintf('correct_block_%d_to_timeline_%d.npy', ...
    db(k).cwExp,db(k).tlExp)]);
aln.pas2tl = readNPY([alignFolder sprintf('correct_block_%d_to_timeline_%d.npy', ...
    db(k).passiveExp,db(k).tlExp)]); % these are in 'timeline' time
aln.noise2tl = readNPY([alignFolder sprintf('correct_block_%d_to_timeline_%d.npy', ...
    db(k).noiseExp,db(k).tlExp)]);

%% load neural data
chanMap = readNPY([ksFolder 'channel_map.npy']);
spks = loadKSdir(ksFolder); % loads all sorts of things (check out the function)

filenamestruct = dir([dataFolder spks.dat_path]);
dataTypeNBytes = numel(typecast(cast(0, spks.dtype), 'uint8')); % determine number of bytes per sample
spks.nSampDat = filenamestruct.bytes/(spks.n_channels_dat*dataTypeNBytes); % to get number of samples

spks.dudChans = find(~ismember(0:384,chanMap)); % indices of the reference channels
spks.liveChans = find(ismember(0:384,chanMap));

[spks.spikeAmps, spks.spikeDepths, spks.templateDepths, ~, spks.tempsUnW, ~, ~] = ...
    templatePositionsAmplitudes(spks.temps, spks.winv, spks.ycoords,  ...
    spks.spikeTemplates, spks.tempScalingAmps);

if any(spks.cgs==2)
    spks.singles = ismember(spks.clu,spks.cids(spks.cgs==2 | spks.cgs==1));
    spks.isGood = ismember(spks.clu,spks.cids(spks.cgs==2));
else
    spks.singles = ismember(spks.clu,spks.cids);
    spks.isGood = [];
end
