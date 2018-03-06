function [ksFolders, dataDir, alnDir, blkDir, alfDir] = expDirs(mouse,date,thisTag,whichServer)
% [ksFolders, dataDir, alnDir, blkDir, alfDir] = expDirs(mouse,date,thisTag,whichServer)
% 
% allows for kilosort files to be located in a different place than data

p = dat.paths();

expFolders = [mouse '\' date '\'];

ksFolders = [expFolders '\ephys_' thisTag '\sorting\'];

dataRoot = [p.([whichServer 'Repository']) '\'];
dataDir = [dataRoot expFolders 'ephys_' thisTag '\'];
alnDir = [dataRoot expFolders 'alignments\'];

if strcmp(whichServer,'old')
    blkDir = ['\\ZSERVER\Data\expInfo\' expFolders];
else
    blkDir = ['\\zubjects.cortexlab.net\Subjects\' expFolders];
end

alfDir = [dataRoot expFolders 'alf\'];