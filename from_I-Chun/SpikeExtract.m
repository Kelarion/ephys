% Spikes = SpikeExtract(DatFileName, SpikeTimes, NumChannels,
%           Channels2Extract, BeforeSamps, AfterSamps)
%
%  will produce an array of extracted waveforms from a .dat file
%  at times specified by SpikeTimes
%
%  First array dimension : channel number
%  Second dimension: time within spike waveform
%  Third dimension: spike number
%
%  SpikeTimes gives the spike times in samples - i.e. a .res file
%  it can be an array or a file name.
%
%  Assumes input is centered on 2048.
%
%  SamplingRate should be in Hz
%
% NB the spikes are NOT resampled and realigned on the peak of the intracellular channel
% if you want that, use SpikeExtractRealign
%
% BeforeSamps and AfterSamps give the number of samples to store on either side of the spike
% (default values = 100)

function Spikes = SpikeExtract(DatFileName, SpikeTimesInput, NChannels, ...
					Channels2Extract, BeforeSamps, AfterSamps)


% Parameters: store number of samples before and after
WaveSamps = BeforeSamps + AfterSamps + 1;

% if SpikeTimesInput is a file name, then load it  ... otherwise assume it is
% an array of spike times.
if (isstr(SpikeTimesInput))
	SpikeTimes = load(SpikeTimesInput);
else
	SpikeTimes = SpikeTimesInput;
end;

SpikeTimes = SpikeTimes(:);
nSpikes = size(SpikeTimes, 1);

fp = fopen(DatFileName, 'r');

% set up array
nChannels2Extract = length(Channels2Extract);
Spikes = zeros(nChannels2Extract, WaveSamps, nSpikes);

% load up data

fprintf('Loading data...\n');

for i = 1:nSpikes
	SkipThis = 0;
	% go to correct part of wave file

	StartSamp = round(SpikeTimes(i) - BeforeSamps);
	status = fseek(fp, StartSamp*NChannels*2, 'bof');
		
	% if error, stop loading spikes
	if (status~=0)
		fprintf('Warning: could not load spike %d\n', i);
		SkipThis = 1;
%		ferror(fp)
%		break;
	end
		
	% Read in buffer
	WaveData = fread(fp, [NChannels WaveSamps], 'int16');

	if any(size(WaveData) ~= [NChannels, WaveSamps])
		fprintf('Warning: could not load spike %d\n', i);
		SkipThis = 1;
	end
	
	% Extract the channels we want, and lose DC offset
	if ~SkipThis
		Spikes(:, :, i) = WaveData(Channels2Extract,:) - 2048;
	else
		Spikes(:,:,i) = 0;
	end
end
	
