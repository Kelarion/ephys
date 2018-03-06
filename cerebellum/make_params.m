% A script for dealing with the quirks of JRClust, organising everything
% into a useable structure "params_all"
replace_existing = 0; % set to true if you want to overwrite

clear db
janelia_db;
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
        [traj, t_spk, t_frame, param] = process_neuropixels_wheel_MATTEO(db{iSet}.trk_files, ...
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
