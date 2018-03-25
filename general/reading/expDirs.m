function [dsetFolders, dataDir, alnDir, blkDir, alfDir] = expDirs(mouse,date,tag,whichServer)
% [dsetFolders, dataDir, alnDir, blkDir, alfDir] = expDirs(mouse,date[,thisTag,whichServer])
% 
% returns the default folders and directories for different kinds of data

if nargin<4, whichServer = 'main'; end
if nargin<3, thisTag = ''; 
else, thisTag = ['ephys_' tag]; end

p = dat.paths();

expFolders = [mouse '\' date '\'];
dsetFolders = [expFolders thisTag '\'];

dataRoot = [p.([whichServer 'Repository']) '\'];
dataDir = [dataRoot dsetFolders '\'];
alnDir = [dataRoot expFolders 'alignments\'];

if strcmp(whichServer,'old')
    blkDir = ['\\ZSERVER\Data\expInfo\' expFolders];
else
    blkDir = ['\\zubjects.cortexlab.net\Subjects\' expFolders];
end

alfDir = [dataRoot expFolders 'alf\'];