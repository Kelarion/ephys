function events = thresholdMinDur(signal,thr,dur)
% events = thresholdMinDur(signal,thr,dur)
%
% Finds the times when a signal is above a threshold value for a certain
% minimum duration. Return as a logical vector the same size as signal.

signal = signal(:)';

aboveThr = signal > thr; % when is it above threshold
evs = conv(aboveThr,ones(dur,1),'same') == dur; % surely theres a better way to do this
events = logical(conv(evs,ones(dur,1),'same')); % when the signal was > thr for the entire bin

end