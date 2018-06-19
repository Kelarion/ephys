

function [stimTimeInds, stimPositions, stimArray] = computeSparseNoiseSignals(block)
% stimTimes are INDICES of which frame showed the stimulus in
% stimPositions. So you should index your photodiode times with it.


sw = block.stimWindowUpdateTimes;

ev = block.events;
stimArray = ev.stimuliOnValues>0;
stimArray = reshape(stimArray, size(stimArray,1), size(stimArray,2)/numel(sw), []);

% these are defined in vis.checker4, which you can see is used in 
% \\zserver\data\expInfo\Erlanger\2016-09-05\2\2016-09-05_2_Erlanger_expDef.m
% Am asking Chris about how to define these within my function itself. They
% should be global parameters.
xRange = [-132 132];
xPos = linspace(xRange(1), xRange(2), size(stimArray,2)+1);
xPos = xPos(1:end-1)+mean(diff(xPos))/2;
yRange = [-36 36];
yPos = linspace(yRange(1), yRange(2), size(stimArray,1)+1);
yPos = yPos(1:end-1)+mean(diff(yPos))/2;

stimArrayZeroPad = cat(3,zeros(size(stimArray,1), size(stimArray,2),1), stimArray);
stimTimeInds = {[], []};
stimPositions = {[], []};
for x = 1:size(stimArray,1)
    for y = 1:size(stimArray,2)
        stimEventTimes{x,y,1} = find(stimArrayZeroPad(x,y,1:end-1)==0 & ...
            stimArrayZeroPad(x,y,2:end)==1); % going from grey to white
        stimEventTimes{x,y,2} = find(stimArrayZeroPad(x,y,1:end-1)==0 & ...
            stimArrayZeroPad(x,y,2:end)==-1); % going from grey to black
        stimTimeInds{1} = [stimTimeInds{1}; stimEventTimes{x,y,1}];
        stimTimeInds{2} = [stimTimeInds{2}; stimEventTimes{x,y,2}];
        
        nEv = length(stimEventTimes{x,y,1});
        stimPositions{1} = [stimPositions{1}; yPos(x)*ones(nEv,1) xPos(y)*ones(nEv,1)];
        nEv = length(stimEventTimes{x,y,2});
        stimPositions{2} = [stimPositions{2}; yPos(x)*ones(nEv,1) xPos(y)*ones(nEv,1)];
    end
end

