function spks = loadJRCdir(jrcDir,getInds)
% spks = loadJRCdir(jrcDir[,getInds])
%
% load all neural data and some meta-data from the output of JRClust.
% jrcDir should be the folder with the *.bin file, and should contain a
% subfolder \clustered\ which has the *_jrc.mat file.
%
% Note: setting 'getInds' to true assumes that 'ind.mat' already exists. 

if nargin <2, getInds = 0; end
mFiles = dir([jrcDir '*.meta']);
neurf  = dir([jrcDir '*.bin']);
if length(mFiles) == 1
    metafile = mFiles(1).name; 
    neurofile = neurf(1).name;
elseif contains(mFiles(1).name,'imec')
    metafile = mFiles(contains({mFiles.name},'ap')).name; 
    neurofile = neurf(contains({neurf.name},'ap')).name;
end 
spks.dat_path = neurofile;
meta_text = fileread([jrcDir metafile]);

[~, ind_nchan] = regexp(meta_text,'nSavedChans=\d');
spks.n_channels_dat = str2num(meta_text(ind_nchan+(0:2)));

[~, ind_prb] = regexp(meta_text,'typeThis='); % adding some flexibiltiy for probe type
probe = meta_text(ind_prb + (1:2)); % will return 'im' for imec, 'ni' for nidq

[~, ind_fs] = regexp(meta_text,[probe 'SampRate=\d']); % it's 'niSampRate' or 'imSampRate' in the meta file
spks.sample_rate = str2num(meta_text(ind_fs+(0:4)));
fs_neuro = spks.sample_rate;

if getInds
    tmp = load([jrcDir '\ind']); % assuming this exists
    event_inds = tmp.ind; % made by get_event_inds[_neuropixels].m
    i_pulse = event_inds{1}; % only for Britton's data
    spks.auxChannel = event_inds(2:end);
    spks.i_pulse = i_pulse;
end

jrc_file = ls([jrcDir '\clustered\*_jrc.mat']);
s = load([jrcDir '\clustered\' jrc_file]); % get the neural data
nClu = s.S_clu.nClu;
t_spk = cell(nClu,1);   notes = cell(nClu,1);
amps = cell(nClu,1);    cids = cell(nClu,1);
for j = 1:nClu % load all info
    cluInds = s.S_clu.cviSpk_clu{j}; % indices of cluster j
    t_spk{j} = double(s.viTime_spk(cluInds))*1000/fs_neuro;
    notes{j} = s.S_clu.csNote_clu{j};
    thisCID = mode(s.S_clu.viClu(cluInds));
    cids{j} = ones(s.S_clu.vnSpk_clu(j),1)*double(thisCID);
    amps{j} = double(s.vrAmp_spk(cluInds));
end
tmp_spk = cell2mat(t_spk);
tmp_cids = cell2mat(cids);
tmp_amps = cell2mat(amps);
[st, inds] = sort(tmp_spk);

dataTypeNBytes = numel(typecast(cast(0, s.P.vcDataType), 'uint8'));
nSampDat = neurf.bytes/(spks.n_channels_dat*dataTypeNBytes);

spks.dtype          = s.P.vcDataType;
spks.nSampDat       = nSampDat;
spks.st             = st; % good old spike times 
spks.clu            = tmp_cids(inds); % cluster ID for each spike
spks.cids           = unique(tmp_cids); % all available cluster IDs
spks.mainChannel    = s.S_clu.viSite_clu; % site is index of channel map
spks.notes          = notes; % 'single' or 'multi'
spks.spikeAmplitude = tmp_amps(inds);
spks.clu_xpos       = s.S_clu.vrPosX_clu; % inferred x position of neuron
spks.clu_ypos       = s.S_clu.vrPosY_clu; % inferred y position of neuron
spks.meanWF         = s.S_clu.tmrWav_raw_clu;
spks.chanMap        = s.P.viSite2Chan; % map from 'site' to row of raw data array
spks.xcoords        = s.P.mrSiteXY(:,1); % x and 
spks.ycoords        = s.P.mrSiteXY(:,2); % y position of those sites


