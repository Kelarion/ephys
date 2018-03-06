function [lf] = getLFP(raw, fs, filt_low)
% Takes a raw signal and gives back the LFP, which is just the signal
% bandpass-filtered between 10 and 50 Hz. 

raw = double(raw);

if nargin < 3
%     filt_low = flt(fs,'lowband');
    c = firpmord([15 20 40 45], [0 1 0], [0.01 0.001 0.01], ...
        fs, 'cell');
    if mod(c{1},2) == 1, c{1} = c{1} + 1; end
    filt_low = firpm(c{:});
end

filtDat = fftfilt(filt_low,raw);
D = floor(length(filt_low)/2);
lf = filtDat(1+D:end);

end