function [WT, F, COI] = cwtnarrow(X,Fs,freqs,varargin)
% [WT, F, COI] = cwtnarrow(X,Fs,freqs[,Name,Value])
% 
% Computes the continuous wavelet transform, between a narrower band than
% determined by the energy of the signal. Specifically, we downsample
% the signal to twice the highest frequency in freqs, setting the upper 
% bound on the cwt, and then choose the range of wavelet scales based on 
% the lowest frequency. See the 'cwt' function for more details.
%
% Inputs: 
%  - X [1,nSamp]:   The signal (narrow-band filter for best results)
%  - Fs [1,1]:      Sampling frequency (in Hz)
%  - freqs [1,2]:   Upper and lower bounds of frequency range
%  - varargin:      All other desired inputs to the 'cwt' function are
%   given as name-value pairs, e.g. ('waveform','amor') to use the Morlet 
%   wavelet. Do NOT specify 'NumOctaves', since that's how this function 
%   sets the lower frequency bound. Look up 'cwt' for default values.
%
% Outputs:
%  - WT [nF,nT]:    Coefficients of the continuous wavelet transform
%  - F [nF,1]:      Frequency of each row in WT
%  - COI [1,nT]:    Indicates where to expect edge effects    

if ~isempty(varargin) % parse varargin
    if mod(length(varargin),2) > 0
        error('Name-value error: not enough input arguments')
    end
    theseNames = {varargin{1:2:end}};
    theseVals = {varargin{2:2:end}};
    specialArgs = false(length(varargin),1);
    if any(contains(lower(theseNames),'wavename'))
        wavename = theseVals{contains(lower(theseNames),'wavename')};
        specialArgs = specialArgs | contains(lower(theseNames),'wavename');
    else
        wavename = 'morse';
    end
    if any(contains(lower(theseNames),'voicesperoctave'))
        Nv = theseVals{contains(lower(theseNames),'voicesperoctave')};
        specialArgs = specialArgs | contains(lower(theseNames),'voicesperoctave');
    else
        Nv = 10; % cwt default
    end
    if any(contains(lower(theseNames),'numoctaves'))
        warning('\nHey, what''d I tell you about specifying NumOctaves?\n')
        specialArgs = specialArgs |contains(lower(theseNames),'numoctaves');
    end
    nOtherArgs = length(varargin) - 2*sum(specialArgs);
    otherargs = cell(nOtherArgs,1);
    otherargs{1:2:end} = theseNames{~specialArgs};
    otherargs{2:2:end} = theseVals{~specialArgs};
else
    wavename = 'morse';
    Nv = 10;
    otherargs = {};
end

nyqFs = freqs(1)*2;
dsFactor = round(Fs/nyqFs);
newX = X(1:dsFactor:end);

a_max = nyqFs/freqs(2); % dilation of the lowest-frequency wavelet
No = log2(a_max) + (1/Nv); % based on a_max = 2^(No - 1/Nv)

[WT, F, COI] = cwt(newX,wavename,nyqFs,'NumOctaves',No,otherargs);

end


