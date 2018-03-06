i = 1;
all_raw = {};
all_lief = {};
all_phase = {};
tfrq = 31; % in Hz
thr = 2;
smothwin = 0.25; % how many seconds to smooth over

for iSet = 1:length(db)
    metafile = [db{iSet}.neuro_file(1:(end-3)) 'meta'];
    meta_text = fileread(metafile);
    [~, bar] = regexp(meta_text,'fileTimeSecs=\d');
    [foo, ~] = regexp(meta_text,'firstSample=\d');
    fileSecs = str2num(meta_text(bar:foo-1));
    fileSamp = round(25000*fileSecs);
    
    singleUnits = contains(param_all{iSet}.notes,'single');
    goodClus = all_spks{iSet}(singleUnits,:);
    for iClu = find(singleUnits)
        actual_site = channels(param_all{iSet}.site(iClu)); % because ephys sucks
        
        raw_channel = double(load_channel('*int16',db{iSet}.neuro_file, ...
            param_all{iSet}.nchan, actual_site, 1, fileSamp));
        all_raw{i} = raw_channel;
        
        lief = getLFP(raw_channel,25000,filt_low);
        all_lief{i} = lief;
        angus = angle(hilbert(lief));
 
        [ron roff]  = findFreq(raw_channel,lief, tfrq,thr,smothwin,0);
        phases = [];
        for iRip = 1:length(ron)
            slice = logical(all_spks{iSet}(iClu,ron(iRip):roff(iRip)));
            phases = [phases angus(slice)];
        end
        all_phase{i} = phases;
        
        i = i + 1;
    end
    
end


for i = 1:length(all_raw)
   figure;polarhistogram(all_phase{i},25,'facecolor','r')
   title(sprintf('cell %d', i))
end