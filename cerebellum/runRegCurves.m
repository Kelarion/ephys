%% Parameters
path(goldenpath) % correct the MatLab path
dset2plot           = 3; % change from this to plot a different dataset
col_parts           = [1 .4 0; 1 0 .4; 0 .4 1; .4 0 1; 0 .7 0];
col_parts_light     = [1 .5 .2; 1 .5 .7; .5 .7 1; .7 .5 1; .3 .7 .3];
replace_existing    = 0; % whether to replace existing files

%% Make the necessary files
clear db
make_db;
for iSet = 1:length(db)
    
    neuro_files_all{iSet} = db{iSet}.neuro_file;
    trk_files_all{iSet} = db{iSet}.trk_files;
    calib_files_all{iSet} = db{iSet}.calib_files; 
    
    slashes = strfind(db{iSet}.neuro_file,'\'); % check for ind file
    neurodir = db{iSet}.neuro_file(1:slashes(end));
    indname = fullfile(neurodir,'ind.mat');
    if ~exist(indname,'file') || replace_existing
        fprintf('Making ind for dataset %d \n', iSet)
        if contains(db{iSet}.neuro_file, 'imec') % depends on the nomenclature of the file! should contain 'imec'
            ind = get_event_ind_neuropixels(db{iSet}.neuro_file, []);
        else
            ind = get_event_ind(db{iSet}.neuro_file,66,65:66,[],0); % make and save ind file
        end
        save(indname,'ind');
    end
    
    clusdir = [neurodir 'clustered'];
    trajname = fullfile(clusdir,'traj.mat');
    spkname = fullfile(clusdir,'t_spk.mat');
    framename = fullfile(clusdir,'t_frame.mat');
    paramname = fullfile(clusdir,'param.mat');
    
    if ~exist(trajname,'file') || replace_existing
        fprintf('processing dataset %d \n',iSet)
        [traj t_spk t_frame param] = process_neuropixels_wheel_MATTEO(db{iSet}.trk_files, ...
                                                    db{iSet}.calib_files,db{iSet}.neuro_file);
        
        save(trajname,'traj');
        save(spkname,'t_spk');
        save(framename,'t_frame');
        save(paramname,'param');
        
        traj_all{iSet} = traj;
        t_spk_all{iSet} = t_spk;
        t_frame_all{iSet} = t_frame;
        param_all{iSet} = param;
    else
        foo = load(trajname);
        traj_all{iSet} = foo.traj;
        foo = load(spkname);
        t_spk_all{iSet} = foo.t_spk;
        foo = load(framename);
        t_frame_all{iSet} = foo.t_frame;
        foo = load(paramname);
        param_all{iSet} = foo.param;
    end

end

%% Run Britton's code

disp('Running regression')
runBrittonRegression;

%% make ETAs

for iSet = 1:length(db)
    if ~contains(db{iSet}.neuro_file, 'imec') % for now, only nidq data has opto
        
        bslash = strfind(param_all{iSet}.rawdir,'\'); % need to do this for now
        uscore = strfind(param_all{iSet}.rawdir,'_'); % make 'make_db' more sensible later
        name = param_all{iSet}.rawdir(bslash(end-1)+1:bslash(end)-1); 
        date = param_all{iSet}.rawdir(bslash(end)+1:uscore(end)-1);
        depth = param_all{iSet}.rawdir(uscore(end)+1:end);
        
        metafile = [db{iSet}.neuro_file(1:(end-3)) 'meta'];
        meta_text = fileread(metafile);
        [~, bar] = regexp(meta_text,'fileTimeSecs=\d');
        [foo, ~] = regexp(meta_text,'firstSample=\d');
        fileSecs = str2num(meta_text(bar:foo-1));
        fileSamp = round(25000*fileSecs);
        
        load([param_all{iSet}.rawdir '\ind']);
        prb_file = dir([param_all{iSet}.rawdir '\*.prb']);
        probe_specs = fileread([prb_file.folder '\' prb_file.name]);
        upto = strfind(probe_specs, 'shank');
        eval(probe_specs(1:upto(2)-1)); % forced to do this because of weird file types
        
        istart = t_frame_all{iSet}(1)*25; % only load things when we were collecting video
        nsamp = t_frame_all{iSet}(end)*25 - istart;
        
%         raw_aux = load_channel('*int16',db{iSet}.neuro_file, ... 
%                                param_all{iSet}.nchan, 66, 1, Inf); % get opto times
%         opto = cleanAux(raw_aux);
%         dopto = diff(opto); % other useful things to know
%         opton = find(dopto > 0);
%         optoff = find(dopto < 0);
%         volleys = zeros(size(opto)); % convert opto pulses to volleys
%         end_volley = diff(find(dopto > 0)) < 10^4;
%         for iOp = 1:length(opton)-1
%             if end_volley(iOp), volleys(opton(iOp):optoff(iOp+1)) = 1; end
%         end
%         voloff = diff([volleys 0]) < 0; % end of volley
        
        % plot event-triggered averages
%         fprintf('plotting dataset %d \n',iSet)
%         fig = figure('Position', [49 49 1822 1068], 'Visible', 'off');
        nClu = length(t_spk_all{iSet});
        dims = factor(nClu); % how best to organise subplots
        h = ceil(length(dims)/2);
        nrow = prod(dims(1:h));
        
        winder = 1.5*param_all{iSet}.fs;
        all_spks{iSet} = zeros(nClu,fileSamp);
        for iClu = 1:nClu
            spks_samp = round(t_spk_all{iSet}{iClu}.*25); % convert to samples
            all_spks{iSet}(iClu,spks_samp) = param_all{iSet}.site(iClu);
        end
    end
end

%% Opto ETAs
clus2plot = [7 10 11 18];
plotDir = 'E:\matteodata\poster\';
for i = 1:length(clus2plot)
    iClu = clus2plot(i);
    spks_samp = round(t_spk_all{iSet}{iClu}.*25); % convert to samples
    all_spks{iSet}(iClu,spks_samp) = param_all{iSet}.site(iClu);
    
    gg = gausswin(0.2*param_all{iSet}.fs); % gaussian window
    %             gg = gg./sum(gg); % normalised to unit integral
    smoth = conv(logical(all_spks{iSet}(iClu,:)),gg,'same'); % gives firing rate
    
    actual_site = channels(param_all{iSet}.site(iClu)); % because ephys sucks
    
    raw_channel = double(load_channel('*int16',db{iSet}.neuro_file, ...
                                   param_all{iSet}.nchan, actual_site, 1));
    
    [eta, peris] = etrigav([smoth; all_spks{iSet}(iClu,:); opto],voloff, ...
        1.5*param_all{iSet}.fs, 0.75*param_all{iSet}.fs);
    
    t = linspace(-1, 1, 1.5*param_all{iSet}.fs);
    spks_peri = cell(1,sum(voloff)); % spike times around opto
    for iev = 1:sum(voloff)
        spks = logical(peris(iev,:,2));
        spks_peri{iev} = t(spks);
    end
    
    op = [diff(peris(2,:,3)) 0]; % opto times
    opTimes = t(logical(op));
    
    %             ax = subaxis(nrow, nClu/nrow, iClu,'Spacing',0.03, 'Padding', 0, 'Margin', 0.03);
    fig = figure('Position' ,[680 303 888 795]); ax = axes;
    title(ax, sprintf('Cell %d on site %d', iClu, param_all{iSet}.site(iClu)))
    
    rasterdot(spks_peri,'axes',ax,'tick_height',0.5,'color',[0.2 0.2 0.2],'spacing',0.8);
    hold on; ylim auto; % plot the spikes and the smoothed PSTH
    plot(t,eta(1,:) + sum(voloff)*0.7,'linewidth',3,'color','k')
    
    opFill = reshape([opTimes; opTimes],1,length(opTimes)*2); % plot opto times
    opHeight = repmat([ax.YLim fliplr(ax.YLim)], 1, length(opFill)/4);
    fill(opFill, opHeight, [0.1 0.4 1], 'EdgeAlpha',0); alpha 0.3
    axis(ax,'off')
    
    hgexport(fig,[plotDir 'Optotag_cell' num2str(iClu)],...
        hgexport('factorystyle'), 'Format', 'epsc');
end

%% Spike-triggered averages

for iSet = 1:length(db)
    bslash = strfind(param_all{iSet}.rawdir,'\'); % need to do this for now
    uscore = strfind(param_all{iSet}.rawdir,'_'); % make 'make_db' more sensible later
    name = param_all{iSet}.rawdir(bslash(end-1)+1:bslash(end)-1);
    date = param_all{iSet}.rawdir(bslash(end)+1:uscore(end)-1);
    depth = param_all{iSet}.rawdir(uscore(end)+1:end);
    
    load([param_all{iSet}.rawdir '\ind']);
    prb_file = dir([param_all{iSet}.rawdir '\*.prb']);
    probe_specs = fileread([prb_file.folder '\' prb_file.name]);
    upto = strfind(probe_specs, 'shank');
    eval(probe_specs(1:upto(2)-1)); % forced to do this because of weird file types
    
    fprintf('plotting dataset %d \n',iSet)
    fig = figure('Position', [49 49 1822 1068], 'Visible', 'off');
    nClu = length(t_spk_all{iSet});
    dims = factor(nClu); % how best to organise subplots
    h = ceil(length(dims)/2);
    nrow = prod(dims(1:h));
    
    dispos = zeros(nClu, 5); % will collect total displacement of each limb
    peakvel = zeros(nClu, 5); % the peak spike-triggered speed for each limb
    peaklag = zeros(nClu, 5); % the tiem of peak from 0
    for iClu = 1:nClu
        [sta peris] = sta_velocity(traj_all{iSet},t_spk_all{iSet}{iClu}', ... 
                                    t_frame_all{iSet}, [-2000 2000], 1);
        
        
        
%         % for right forelimb
%         spead = sqrt(diag(sta.sta_v_RFL_SS'*sta.sta_v_RFL_SS)); % convert to speed (sqrt(<v(t),v(t)>)
%         dispos(iClu,1) = sum(spead); % and get the total spike-triggered displacement 
%         peakvel(iClu,1) = max(spead);
%         peaklag(iClu,1) = find(spead == max(spead)) - floor(length(spead)/2);
%         
%         % for left forelimb
%         spead = sqrt(diag(sta.sta_v_LFL_SS'*sta.sta_v_LFL_SS)); % convert to speed (sqrt(<v(t),v(t)>)
%         dispos(iClu,2) = sum(spead); % and get the total spike-triggered displacement 
%         peakvel(iClu,2) = max(spead);
%         peaklag(iClu,2) = find(spead == max(spead)) - floor(length(spead)/2);
%         
%         % for right hindlimb
%         spead = sqrt(diag(sta.sta_v_RHL_SS'*sta.sta_v_RHL_SS)); % convert to speed (sqrt(<v(t),v(t)>)
%         dispos(iClu,3) = sum(spead); % and get the total spike-triggered displacement 
%         peakvel(iClu,3) = max(spead);
%         peaklag(iClu,3) = find(spead == max(spead)) - floor(length(spead)/2);
%         
%         % for left hindlimb
%         spead = sqrt(diag(sta.sta_v_LHL_SS'*sta.sta_v_LHL_SS)); % convert to speed (sqrt(<v(t),v(t)>)
%         dispos(iClu,4) = sum(spead); % and get the total spike-triggered displacement 
%         peakvel(iClu,4) = max(spead);
%         peaklag(iClu,4) = find(spead == max(spead)) - floor(length(spead)/2);
%         
%         % for the butt
%         spead = sqrt(diag(sta.sta_v_Bt_spk'*sta.sta_v_Bt_spk)); % convert to speed (sqrt(<v(t),v(t)>)
%         dispos(iClu,5) = sum(spead); % and get the total spike-triggered displacement 
%         peakvel(iClu,5) = max(spead);
%         peaklag(iClu,5) = find(spead == max(spead)) - floor(length(spead)/2);
%         
    end
    
    
end

%% 

clus2plot = [7 10 11 18];
plotDir = 'E:\matteodata\poster\';
for i = 1:length(clus2plot)
    iClu = clus2plot(i);
    
    
    
    hgexport(fig,[plotDir 'Optotag_cell' num2str(iClu)],...
        hgexport('factorystyle'), 'Format', 'epsc');
end
