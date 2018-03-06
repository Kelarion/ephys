function aln = loadAlignments(alignFolder,thisTag,tlExp,varargin)
% aln = loadAlignments(alignFolder,thisTag,tlExp[,'expName',expNum])
%
% output will contain alignment between given ephys tag and the reference
% probe. subsequent inputs load the alignment between the given exp and the
% reference time (e.g. timeline to reference ephys, or block_3 to timeline)

if ~isempty(varargin)
    if mod(length(varargin),2) > 0
        error('Name-value error: Need an ''expNum'' for each ''expName''')
    end
    theseExps = varargin(1:2:end);
    expNums = varargin(2:2:end);
else
    theseExps = {};
end

alignToRef = dir([alignFolder 'correct_timeline_*.npy']);
[~, refInd] = regexp(alignToRef.name,'ephys_');
aln.refTag = alignToRef.name(refInd+[1:2]);

aln.thisTag = thisTag;
if strcmp(thisTag,aln.refTag)
    aln.tag2ref = [1 0];
else
    aln.tag2ref = readNPY([alignFolder sprintf('correct_ephys_%s_to_ephys_%s.npy', ...
        thisTag,aln.refTag)]); % correct from ref ephys to current ephys
end
aln.tl2ref = readNPY([alignToRef.folder '\' alignToRef.name]);

for iExp = 1:length(theseExps)
    fn = sprintf('correct_block_%d_to_timeline_%d.npy', expNums{iExp},tlExp);
    if exist([alignFolder fn],'file')
        aln.([theseExps{iExp} '2tl']) = readNPY([alignFolder fn]);
    else
        try
            fn = sprintf('mpep_%d_onsets_in_timeline_%d.npy', expNums{iExp},tlExp);
            alnFunc = [1; readNPY([alignFolder fn])];
            aln.([theseExps{iExp} '2tl']) = alnFunc;
        catch
            continue
        end
    end
end