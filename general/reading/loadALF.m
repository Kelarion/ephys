function beh = loadALF(alfDir,varargin)
% beh = loadALF(alfDir[,extraFields])
%
% default fields are:
% {'cwFeedback' 'cwGoCue' 'cwResponse' ... 
% 'cwStimOn' 'cwTrials' 'licks' 'passiveBeep' 'passiveStimOn'  ... 
% 'passiveTrials' 'passiveValveClick' 'wheel' 'wheelMoves'}

fields = {'cwFeedback' 'cwGoCue' 'cwResponse' ... 
    'cwStimOn' 'cwTrials' 'licks' 'passiveBeep' 'passiveStimOn'  ... 
    'passiveTrials' 'passiveValveClick' 'wheel' 'wheelMoves'};

if ~isempty(varargin)
    fields(end:end+length(varargin)) = varargin;
end
fields = unique(fields);

beh = struct;
for iField = 1:length(fields)
    fils = dir([alfDir filesep fields{iField} '.*.npy']);
    for iFil = 1:length(fils)
        subfield = split(fils(iFil).name,[fields{iField} '.']);
        subfield = subfield{end}(1:regexp(subfield{end},'.npy')-1);
        beh.(fields{iField}).(subfield) = readNPY([alfDir fils(iFil).name]);
    end
end

end