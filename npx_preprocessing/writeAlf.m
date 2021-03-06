% Script for writing structures in the alf folder, taken more or less
% wholesale from Nick's behavToALF script but adapted to my needs.

% by default skip datasets with pre-existing structures
if ~exist('overwrite','var'), overwrite = false; end 
%% sparse noise
for k = 1:length(db)
    
    mouseName = db(k).mouse_name;
    thisDate = db(k).date;
    
    root = getRootDir(mouseName, thisDate);
    
    alignDir = fullfile(root, 'alignments');
    alfDir = getALFdir(mouseName, thisDate);
    
    if ~exist(alfDir,'dir')
        mkdir(alfDir);
    elseif exist(fullfile(alfDir, 'sparseNoise.positions.npy'),'file') && (~overwrite)
        continue
    end
    
    % determine which block is the correct one for mapping data
    
    rootExp = dat.expFilePath(mouseName, thisDate, 1, 'Timeline', 'master');
    expInf = fileparts(fileparts(rootExp));
    
    d = dir(fullfile(expInf, '*'));
    expNums = cell2mat(cellfun(@str2num, {d(3:end).name}, 'uni', false));
    
    rfExpNums = [];
    for e = 1:length(expNums)
        % if block, load block and get stimWindowUpdateTimes
        dBlock = dat.expFilePath(mouseName, thisDate, expNums(e), 'block', 'master');
        if exist(dBlock, 'file')
            load(dBlock);
            if isfield(block, 'expDef') && ~isempty(strfind(block.expDef, 'sparseNoiseAsync_NS'))
                rfExpNums(end+1) = e;
            end
        end
    end
    
    % determine which timeline to use
    % defined as which one we have the sync for
    
    for e = length(rfExpNums):-1:1 % search in reverse order, to choose the last if multiple
        d = dir(fullfile(alignDir, sprintf('block_%d_sw*', rfExpNums(e))));
        if ~isempty(d) % found an alignment between this rf and a timeline
            q = sscanf(d.name, 'block_%d_sw_in_timeline_%d.npy');
            rfExpNum = q(1); tlExpNum = q(2);
            fprintf(1, 'using rfExpNum %d and tlExpNum %d\n', rfExpNum, tlExpNum);
            break
        end
    end
    
    % load sync information
    
    tlToMasterFile = dir(fullfile(alignDir, ...
        sprintf('correct_timeline_%d_to_ephys_*.npy', tlExpNum)));
    bTLtoMaster = readNPY(fullfile(alignDir,tlToMasterFile.name));
    
    % load experiment info
    load(dat.expFilePath(mouseName, thisDate, rfExpNum, 'block', 'master'))
    
    stimArrayTimes = readNPY(fullfile(alignDir, ...
        sprintf('block_%d_sw_in_timeline_%d.npy', rfExpNum, tlExpNum)));
    
    stimArrayTimes = applyCorrection(stimArrayTimes, bTLtoMaster);
    %
    [stimTimeInds, stimPositions, stimArray] = ...
        computeSparseNoiseSignals(block);
    
    if length(block.stimWindowUpdateTimes)==length(stimArrayTimes)+1 && min(stimTimeInds{1})>1
        % this is the weird case I still can't figure out where sometimes you
        % have to drop the first stimWindowUpdateTimes
        stimTimeInds = cellfun(@(x)x-1, stimTimeInds, 'uni', false);
    end
    
    stimTimes = cellfun(@(x)stimArrayTimes(x), stimTimeInds, 'uni', false);
    
    % write to alf
    alf.writeEventseries(alfDir, 'sparseNoise', stimTimes{1}, [], []);
    writeNPY(stimPositions{1}, fullfile(alfDir, 'sparseNoise.positions.npy'));
    
end
