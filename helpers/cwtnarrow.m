function [WT, F, T, COI] = cwtnarrow(X,Fs,freqs,varargin)
% [WT, F, T, COI] = cwtnarrow(X,Fs,freqs[,Name,Value])
% 
% Wrapper for 'cwt' function.
% Computes the continuous wavelet transform in a narrower band than
% determined by the energy of the signal. Specifically, we downsample
% the signal to twice the highest frequency in freqs, setting the upper 
% bound on the cwt, and then choose the range of wavelet scales based on 
% the lowest frequency. See the 'cwt' function for more details.
%
% Inputs: 
%  - X (1,nSamp):   The signal (narrow-band filter for best results)
%  - Fs (1,1):      Sampling frequency (in Hz)
%  - freqs (1,2):   Upper and lower bounds of frequency range
%  - Optional:      All other desired inputs to the 'cwt' function are
%   given as name-value pairs, e.g. ('waveform','amor') to use the Morlet 
%   wavelet. Do NOT specify 'NumOctaves', since that's how this function 
%   sets the lower frequency bound. Look up 'cwt' for default values.
%
% Outputs:
%  - WT (nF,nT):    Coefficients of the continuous wavelet transform
%  - F (nF,1):      Frequency of each row in WT
%  - T (1,nT):      Time (in sec) of each column
%  - COI (1,nT):    Indicates where to expect edge effects    

if ~isempty(varargin) % parse varargin
    if mod(length(varargin),2) > 0
        error('Name-value error: not enough input arguments')
    end
    nms = varargin(1:2:end); % names 
    rgs = varargin(2:2:end); % args
    include = false(length(varargin)/2,1);
    if any(contains(lower(nms),'wavename'))
        wavename = rgs{contains(lower(nms),'wavename')};
        include = include | contains(lower(nms),'wavename');
    else
        wavename = 'morse';
    end
    if any(contains(lower(nms),'voicesperoctave'))
        Nv = rgs{contains(lower(nms),'voicesperoctave')};
    else
        Nv = 10; % cwt default
    end
    if any(contains(lower(nms),'numoctaves'))
        warning('\nHey, what''d I tell you about specifying NumOctaves?\n')
        include = include |contains(lower(nms),'numoctaves');
    end
    nExtraArgs = length(varargin) - 2*sum(include);
    otherargs = cell(nExtraArgs,1);
    otherargs(1:2:end) = nms(~include);
    otherargs(2:2:end) = rgs(~include);
else
    wavename = 'morse';
    Nv = 10;
    otherargs = {};
end
t = [1:length(X)]/Fs;
nyqFs = freqs(1)*2 + 10;
dsFactor = round(Fs/nyqFs);
newX = X(1:dsFactor:end);
newt = t(1:dsFactor:end);

a_max = nyqFs/(freqs(2)*2); % scale factor of the lowest-frequency wavelet
No = round(log2(a_max) + (1/Nv)); % based on a_max = 2^(No - 1/Nv)

[WT, F, COI] = cwt(newX,wavename,nyqFs,'NumOctaves',No,otherargs{:});
T = newt;

end

