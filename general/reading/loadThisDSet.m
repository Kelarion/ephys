function varargout = loadThisDSet(db,whichDSet,whatThings)
% thing = loadThisDSet(db,whichDSet,whatThings)
%
%

k = whichDSet(end);
thisTag = db(k).tags{whichDSet(end-1)};

[dsetFolders, dataDir, alnDir, blkDir, alfDir] = ...
    expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
if ~isempty(db(k).ksRoot)
    ksDir = [db(k).ksRoot dsetFolders '\sorting\'];
else
    ksDir = [dataDir '\sorting\'];
end

n = 1;
if any(contains(whatThings,'align'))
    aln = loadAlignments(alnDir,thisTag,db(k).tlExp,'noise',db(k).noiseExp);
    varargout{n} = aln;
    n = n+1;
end
if any(contains(whatThings,'spikes'))
    spks = loadNeuralData(ksDir,dataDir);
    varargout{n} = spks;
    n = n+1;
end
if any(contains(whatThings,'alf'))
    beh = loadALF(alfDir);
    varargout{n} = beh;
    n = n+1;
end
if any(contains(whatThings,'blocks'))
    blk = loadBlocks(blkDir,db(k).tlExp,'cw',db(k).cwExp,'pas',db(k).passiveExp);
    varargout{n} = blk;
end