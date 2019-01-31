function [scTop, scBottom] = wheresCortex(lfpName,algo,algoInfo)
% [scTop, scBottom] = wheresCortex(lfpName,algo[,algoInfo])
% 
% Finds the approximate border of cortex and SC based on the LFP. Assumes
% this is a Neuropixels option 3 or 1 probe.
%
% Inputs: 
%   - lfpName [char]: Full path to LFP *.bin file.
%   - algo [char]: How to find the border. Either 'corrs' or 'snrf':   
%       > 'corrs' computes the correlation matrix and looks for
%       block-diagonal structure, assuming that the closest block to the
%       surface corresponds to cortical channels. This is slow.
%       > 'snrf' uses the mean response to visual noise, taking the most
%       dorsal channel with a strong visual response as the surface of the
%       SC. If you have the sparse noise info, this is much faster.
%   -algoInfo: Depends on algo:
%       > If 'corrs', specify [start, end] times, in seconds
%       > If 'snrf', a cell of {timeCourses, channels} with the time
%       courses and channel indices of all the LFP responses (i.e. first
%       and last output of the snrfLFP function).
%
% Note: this is only guaranteed to work with the specific data pipeline
% used for my SC project in the cortex lab.

maxTimeForCorr = 600; % seconds; keep it < ~10 minutes or the array gets too huge
lfp_thresh = 1.2; % threshold on z-scored LFP response 
dudChans = [37    76   113   152   189   228   265   304   341   380]; %internal refs
q = 1:384;
liveChans = q(~ismember(q,dudChans));

switch algo
    case 'corrs'
        lfpnamestruct = dir(lfpName);
        dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
        
        nSampLFP = lfpnamestruct.bytes/(385*dataTypeNBytes);
        lfpFs = 2500; % assumption for imec
        
        lfpFile = memmapfile(lfpName, 'Format', {'int16', [385 nSampLFP], 'x'});
        
        t0 = round(algoInfo(1)*lfpFs);
        t_max = algoInfo(1)+maxTimeForCorr;
        tfin = (nSampLFP-1)/lfpFs;
        whenCompute = t0:round(min([algoInfo(2) t_max tfin])*lfpFs);
        
        mySlice = double(lfpFile.Data.x(liveChans,whenCompute))'; % start at halfway
        medianSubtractedLFP = mySlice - median(mySlice,2);
        cormat = corrcoef(medianSubtractedLFP);
        
        rows = findDiagBlocks(cormat,3,'svd'); % finds row of the block diagonals
        scTop = max(rows);
    case 'snrf' % best option if you can use it
%         tmp = load(algoInfo); % this requires the snrf struct 
%         snrf = tmp.snrf;      % which is made by getRFs.m
%         timeCourse = snrf.lfp_timecourse;
%         tcChans = snrf.computed_channels;
        timeCourse = algoInfo{1}; % the first output of snrfLFP
        tcChans = algoInfo{2}; % last output of snrfLFP
        [m, ind] = max(-timeCourse,[],2);
        peakTime = ind(m == max(m));
        depthResponse = zscore(timeCourse(:,peakTime));
        respChan = tcChans(thresholdMinDur(-depthResponse,lfp_thresh,2));
        scTop = max(respChan);
        scBottom = min(respChan);
    otherwise
        error('Please give a valid algorithm: either ''corrs'' or ''snrf''')
end

end