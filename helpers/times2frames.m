function lgclIndex = times2frames(eventTimes,frameTimes,enumerate)
% align = times2frames(eventTimes,frameTimes[,enumerate])
% 
% Take a series of event times and convert them to quasi-logical vectors,
% where each entry is non-zero if the event is occuring. 
%
% -eventTimes: vector with timestamps of each event (length nTimes)
% -frameTimes: vector with the time of each sample (length nFrames)
% -enumerate (optional): if 1, then each nonzero entry corresponds to the
%  number of that event.
%
% If eventTimes is a nTimes x 2 matrix, then it is assumed that the first
% row corresponds to the onset, and the second the offset. 

if ~exist('enumerate','var'), enumerate = 0; end
nTimes = size(eventTimes,1);
lgclIndex = zeros(1,length(frameTimes));
isCont = min(size(eventTimes)) - 1; % do we report the event continuously?
if size(eventTimes,2) ~= min(size(eventTimes))
   eventTimes = eventTimes'; % make sure that we index properly
end

for iEv = 1:nTimes
    iStart = find(frameTimes >= eventTimes(iEv,1),1,'first');
    if isCont
        iFin = find(frameTimes >= eventTimes(iEv,2),1,'first');
        if enumerate
            lgclIndex(iStart:iFin) = iEv;
        else
            lgclIndex(iStart:iFin) = 1;
        end
    else
        if enumerate
            lgclIndex(iStart) = iEv;
        else
            lgclIndex(iStart) = 1;
        end
    end
end

end