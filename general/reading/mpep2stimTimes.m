function [stimTimes, stimPositions] = mpep2stimTimes(mpepOnsets,pars,hwInfo,tl2ref)
% [stimTimes, stimPositions] = mpep2stimTimes(mpepOnsets,params,hwInfo,tl2ref)
%
% Extracts the stimulus information for sparse visual noise experiments run
% with mpep. Outputs are arrays organised by stimulus type:
%   {allStims, whiteStims, blackStims}
%
% adapted from Sylvia Schroeder, UCL Cortex Lab

stimFile = str2func(strtok(pars.xfile, '.'));
hwInfo.windowPtr = NaN;

noiseOn = [mpepOnsets(:) ones(size(mpepOnsets(:)))]*tl2ref;

SS = stimFile(hwInfo, pars.pars);
stimFrames = cat(3, SS.ImageTextures{:});

framesPerImage = pars.pars(6,1);
frameTimes = (0 : size(stimFrames, 3)-1) * framesPerImage / hwInfo.FrameRate;
ft = (frameTimes + noiseOn)';

rows = size(stimFrames,1);
cols = size(stimFrames,2);

% white squares
ind = find(stimFrames == 1);
t = ceil(ind / (rows * cols));
ind = mod(ind, (rows * cols));
ind(ind == 0) = rows * cols;
x_wh = ceil(ind / rows);
y_wh = mod(ind, rows);
y_wh(y_wh == 0) = rows;
time_wh = ft(t,:);
% black squares
ind = find(stimFrames == -1);
t = ceil(ind / (rows * cols));
ind = mod(ind, (rows * cols));
ind(ind == 0) = rows * cols;
x_bl = ceil(ind / rows);
y_bl = mod(ind, rows);
y_bl(y_bl == 0) = rows;
time_bl = ft(t,:);

stimPositions = cell(1,3);
stimTimes = cell(1,3);
stimPositions{2} = [repmat(y_wh, size(ft,2), 1) repmat(x_wh, size(ft,2), 1)];
stimTimes{2} = reshape(time_wh, [], 1);
stimPositions{3} = [repmat(y_bl, size(ft,2), 1) repmat(x_bl, size(ft,2), 1)];
stimTimes{3} = reshape(time_bl, [], 1);
stimTimes{1} = [stimTimes{2}; stimTimes{3}];
[stimTimes{1}, order] = sort(stimTimes{1}, 'ascend');
stimPositions{1} = [stimPositions{2}; stimPositions{3}];
stimPositions{1} = stimPositions{1}(order,:);



end
