function blks = loadBlocks(blockDir,tlExp,varargin)
% blks = loadBlocks(blockDir,tlInfo,tlExp[,'expName',expNum])
%
% output will contain timeline, as "blks.tl", and then the block structures
% of all subsequent experiments listed.

if ~isempty(varargin)
    if mod(length(varargin),2) > 0
        error('Name-value error: Need an ''expNum'' for each ''expName''')
    end
    theseExps = varargin(1:2:end);
    expNums = varargin(2:2:end);
else
    theseExps = {};
end

tldir = dir([blockDir '\' num2str(tlExp) '\'  '*Timeline.mat']);
tl = load([tldir.folder '\' tldir.name]);
blks.tl = tl.Timeline;

for iExp = 1:length(theseExps)
    thisDir = dir([blockDir '\' num2str(expNums{iExp}) '\*Block.mat']);
    if isempty(thisDir)
        warning(['\nNo Block found for experiment: ' num2str(expNums{iExp}) ', skipping'])
        continue
    end
    tempblock = load([thisDir.folder '\' thisDir.name]);
    blks.(theseExps{iExp}) = tempblock.block;
end