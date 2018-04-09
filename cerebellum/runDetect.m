%% load
parentDir = [ksRoot '\' db(k).name '\'];
subDir = dir([parentDir '*' db(k).date]);

if length(subDir) > 1
    subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
    if isempty(subDir)
        disp('can''t find it ¯\_(?)_/¯')
    end
end
dataDir = [parentDir subDir.name '\'];

sp = loadJRCdir(dataDir);
rawDat = memmapfile([dataDir sp.dat_path], 'Format',  ...
    {sp.dtype, [sp.n_channels_dat sp.nSampDat], 'x'});





