function bin2hdf(fname,fromHere,toHere, compress)
% Converts .bin files to hdf5 format, to make file I/O much more efficient
% and portable. http://geology.beer/2015/02/10/hdf-for-large-arrays/
%
% To avoid clutter, you can supply the directories of the original *.bin
% (fromHere) and the new hdf5 (toHere).  
%
% The major advantages are chunking of the data (for faster I/O and less 
% RAM usage) and the ease of indexing (no need to work in byte offsets, and 
% all formatting is the same across machines and platforms). 
% 
% The trade-off is that the files can be much larger than binary files. 
% Setting 'compress' to 'on' will aleviate this, but writing the files
% takes an order of magnitude more time (gets worse for larger the data).
%
% Matteo Alleman, 2017. JFRC

if nargin < 3
    toHere = fromHere;
end

if nargin < 4
    compress = 'off';
end

switch compress
    case 'on'
        comp = 1;
    otherwise
        comp = 0;
end

tic
binplace = fullfile(fromHere,fname);
hname = strcat(fname(1:end-3),'h5');
hdfplace = fullfile(toHere,hname);

fprintf('%f : Loading bin file ... \n', toc)
binid = fopen(binplace,'r');
dat = fread(binid,[66 Inf],'*int16');
dims = size(dat);

% fprintf('%f : Creating HDF file ... \n', toc)
h5create(hdfplace,'/data',dims,'ChunkSize',[66 25000],'Deflate',comp)
fprintf('%f : Writing to HDF file ... \n', toc)
h5write(hdfplace,'/data',dat)

fclose(binid);
clear dat
fprintf('%f : Done! \n', toc)

end