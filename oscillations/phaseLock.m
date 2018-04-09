function shiftedEventTimes = phaseLock(eventTimes,oscillation,whichPhase,Fs)
% shiftedEventTimes = phaseLock(eventTimes,oscillation,whichPhase,Fs)
% 
% whichPhase is between [-pi pi].

win = round(Fs*[-0.1 0.1]);

evntInds = round(eventTimes*Fs);
dtheta = abs(angle(hilbert(oscillation)) - whichPhase);

shiftedEventTimes = zeros(length(eventTimes),1);
for iEv = 1:length(eventTimes),ii = evntInds(iEv);
    samps = ii+win(1):ii+win(2);
    [~,dt] = min(dtheta(samps(samps>0)));
    dt = dt + win(1);
    shiftedEventTimes(iEv) = eventTimes(iEv) + ((dt(1))/Fs);
end

