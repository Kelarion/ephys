localRoot = 'C:\DATA\Spikes\';

modelStrings = {'gaussian','DoG','DoUG'};% 'Gabor'};
modelFuncs = {@D2GaussFunctionRot, @D2_DoG_Rot, @D2_DoG_unequal_Rot};% @D2GaborRot};
ipsiFits = {'gaussian'};% 'same'};

params.nFold = 5;
params.nStoch = 20;
params.method = 'direct'; % 'itskov' 

rfparams = struct('useSVD',true,'makePlots',false,'fit2Dgauss',false,'nShuffle',200);

min_for_ipsi = 35;
z_thresh = 2; % threshold on peak RF Z-score (not very useful for me, set to very small value)
nspk_thresh = 720; % minumum number of spikes for analysis (I say they should have at least
                      % been able to spike twice for each stimulus)
saveName = 'all_RF_fits.mat';

% parpool('SpmdEnabled',false); % let the loop keep running if a worker dies

clear db
ephys_RF_db
nDsets = length(db);
%% run 
clu_fits = struct;
wc = cell(nDsets,1); % these are megacells, each dataset is stored as an entry
crf = cell(nDsets,1);
bf = cell(nDsets,1);
bp = cell(nDsets,1);
rs = cell(nDsets,1); % called r-squared because I used only to use that
rfs = cell(nDsets,1);
rfrs =cell(nDsets,1);
null_rs = cell(nDsets,1);
emprcl_rs = cell(nDsets,1);
best_rs = cell(nDsets,1);
clinf = cell(nDsets,1);
dispstat('','init');
parfor k = 1:nDsets
    whichCell = []; % we fill these up with the information from all probes in a dataset;
    clu_rf = {};  % since we don't know the total number of cells across all probes 
    best_fit = {}; % we've got to leave this empty for now
    best_params = {}; % (since I have a threshold # spikes I can't know nClu ahead of time)
    r_squared = []; 
    rf_stats = struct('timeBins',{},'timeCourse',{},'peakZscore',{},'scalar',{},'bsl',{},'z_peak',{});
    rf_rsq = []; % cross-validated empirical model
    flat_rsq = []; % null model, should always be bad
    empirical_rsq = []; % empirical model (un-CV)
    bestmodel_rsq = []; % best model fit on all data (un-CV)
    cluster_info = [];
    for t = 1:length(db(k).tags)
        thisTag = db(k).tags{t};
        
        thisExp =  [db(k).mouse_name '_' db(k).date '_' thisTag];
%         if verbose, dispstat(['Starting ' thisExp],'keepthis','timestamp'); end
        
        % load everything for that probe
        [dsetFolders,dataDir,~,alfDir] = expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [dataDir '\sorting\'];
        
        snname = [localRoot dsetFolders '\sparse_noise_RFs.mat'];
        if exist(snname,'file')
            snrf = loadVar(snname,'snrf');
        else
            disp([thisExp '      doesn''t have an RF file, nevermind'])
            continue
        end
        
        useThese = [snrf.neur_rfstats(:).peakZscore] > z_thresh;
        
        spks = loadNeuralData(ksDir);
        noiseStart = min(snrf.stimTimes_local);
        noiseEnd = max(snrf.stimTimes_local);
        
        bfname =[alfDir thisTag '\borders_' thisTag '.tsv'];
        [scsurf,stratopt] = bord(bfname,snname,spks); % borders from alf or lfp
        
        cluDepth = clusterAverage(spks.clu(spks.isNeuron),spks.spikeDepths(spks.isNeuron));
        cluDepth = cluDepth(:);
        cluAmps = clusterAverage(spks.clu(spks.isNeuron),spks.spikeAmps(spks.isNeuron))*spks.gain;
        cluAmps = cluAmps(:);

        cluNspk = zeros(length(snrf.neur_ID),1);
        for iClu = 1:length(snrf.neur_ID)
            cluNspk(iClu) = sum((spks.clu == snrf.neur_ID(iClu)) & (spks.st>noiseStart & spks.st<noiseEnd));
        end
        
        useThese = useThese(:) & (cluNspk>nspk_thresh);
        visClu = snrf.neur_ID(useThese);
        
        tmp_info = [cluDepth(useThese), ...
            cluDepth(useThese)-scsurf, ...
            cluDepth(useThese)-stratopt, ...
            cluAmps(useThese), ...
            cluNspk(useThese)];
        cluster_info = [cluster_info; tmp_info];
        
        yPos = snrf.YPos; nX = length(yPos);
        xPos = snrf.XPos; nY = length(xPos);
        nMod = length(modelStrings);
        nClu = length(visClu);
        %%
        tmp_rf = cell(nClu,1); % we fill these up looping over clusters
        tmp_rfrsq = zeros(nClu,1);
        tmp_flat = zeros(nClu,1);
        tmp_rsq = nan(nClu,nMod+length(ipsiFits),params.nFold);
        tmp_cids = zeros(nClu,3);
        tmp_sts = struct('timeBins',{},'timeCourse',{},'peakZscore',{},'scalar',{},'bsl',{},'z_peak',{});
        tmp_bf = cell(nClu,1);
        tmp_bfprms = cell(nClu,1);
        tmp_emprcl = zeros(nClu,1);
        tmp_fullfit = zeros(nClu,2);
        dispstat(['Starting ' thisExp ', running CV for ' num2str(nClu) ' cells...'],'keepprev','keepthis','timestamp')
        for iClu = 1:nClu
            CID = visClu(iClu);
%             RF = snrf.neur_rfmap{snrf.neur_ID == CID};
%             sts = snrf.neur_rfstats(snrf.neur_ID == CID);
            
            deezspks = (spks.clu == CID) & (spks.st>noiseStart & spks.st<noiseEnd);
            
            [RF,sts] = sparseNoiseRF_MA(spks.st(deezspks), snrf.stimTimes_local,snrf.stimPosition,rfparams);
            
            rsq = nan(1,nMod,params.nFold);
            mdl = struct('String',{},'func',{},'fitIpsi',{},'ipsitype',{});  
            mnX0 = zeros(1,nMod);
            mnWx = zeros(1,nMod);
            for iMod = 1:nMod
                mdl(iMod).String = modelStrings{iMod};
                mdl(iMod).func = modelFuncs{iMod};
                mdl(iMod).fitIpsi = false;
                mdl(iMod).ipsitype = 'no';
                
                try
                    [rsq(1,iMod,:),x] = cvalGausModel(mdl(iMod),spks.st(deezspks),snrf.stimTimes_local,snrf.stimPosition,params);
                catch thepoxanddieyoubitch
                    disp(thepoxanddieyoubitch.message)
                    error(['this motherfucking buggy piece of shit function fucked up and couldn''t even work on cell ' num2str(CID) ' in ' thisExp])
                end
                mnX0(iMod) = median(round(x(2,:)));
                mnWx(iMod) = median(round(x(3,:)));
                
            end
            
            strawmdl = struct('String','rf','func',[],'fitIpsi',[],'ipsitype',[]);
            score_rf = nanmean(cvalGausModel(strawmdl,spks.st(deezspks),snrf.stimTimes_local,snrf.stimPosition,params));
            
            strawmdl = struct('String','flat','func',[],'fitIpsi',[],'ipsitype',[]);
            score_flat = nanmean(cvalGausModel(strawmdl,spks.st(deezspks),snrf.stimTimes_local,snrf.stimPosition,params));
            
            score_emprcl = empiricalFit(RF,sts,snrf,spks.st(deezspks),params.method);
            
            [~,best_model_for_now] = max(mean(rsq,3));
            % don't bother fitting ipsilateral field if the best fit was too medial
            if abs(mnX0(best_model_for_now)) >= min_for_ipsi
                for jj = 1:length(ipsiFits)
                    mdl(nMod+jj).String = modelStrings{best_model_for_now};
                    mdl(nMod+jj).func = modelFuncs{best_model_for_now};
                    mdl(nMod+jj).fitIpsi = true;
                    mdl(nMod+jj).ipsitype = ipsiFits{jj};
                    
                    rsq(1,iMod+jj,:) = cvalGausModel(mdl(nMod+jj),spks.st(deezspks),snrf.stimTimes_local,snrf.stimPosition,params);
                    
%                     dispstat(['Done with model: ' modelStrings{best_model_for_now} ' with ' ipsiFits{jj} ' ipsilaterally'],'timestamp')
                end
            else
                rsq(1,iMod+1:iMod+length(ipsiFits),:) = NaN;
            end
            
            % get the best fits
            [naan, best_ipsi] = max(mean(rsq(:,nMod+1:end,:),3));
            if isnan(naan)
                best_ipsi = [];
            else
                best_ipsi = best_ipsi+nMod;
            end
            wantedModels = unique([best_model_for_now, best_ipsi]);
            
            [all_fits, prms] = fitTheseModRF(mdl(wantedModels),RF,xPos,yPos,1.5,params.nStoch);
            score_fullmodel = [NaN NaN];
            for ii = 1:size(all_fits,3)
                score_fullmodel(ii) = empiricalFit(all_fits(:,:,ii),sts,snrf,spks.st(deezspks),params.method);
            end
            
            %disp('after fitTheseModRF')
            tmp_rf{iClu} = RF;
            tmp_rsq(iClu,:,:) = rsq;
            tmp_cids(iClu,:) = [CID t k]; %disp('before making tmp_sts')
            tmp_sts(iClu) = sts;% disp('after making tmp_sts')% for debugging
            tmp_bf{iClu} = all_fits;
            tmp_bfprms{iClu} = prms;
            tmp_rfrsq(iClu) = score_rf;
            tmp_flat(iClu) = score_flat;
            tmp_emprcl(iClu) = score_emprcl;
            tmp_fullfit(iClu,:) = score_fullmodel;
        end
        whichCell = [whichCell; tmp_cids];
        clu_rf = cat(1,clu_rf,tmp_rf);
        best_fit = cat(1,best_fit,tmp_bf);
        best_params = cat(1,best_params,tmp_bfprms);
        r_squared = cat(1,r_squared,tmp_rsq);
        rf_stats = cat(1,rf_stats(:),tmp_sts(:));
        rf_rsq = [rf_rsq(:); tmp_rfrsq(:)];
        flat_rsq = [flat_rsq(:); tmp_flat(:)];
        empirical_rsq = [empirical_rsq(:); tmp_emprcl(:)];
        bestmodel_rsq = [bestmodel_rsq; tmp_fullfit];
%         disp('before saving') % for debugging
        fname = [db(k).mouse_name '_' db(k).date '_temp_vars.mat'];
        parsave(fname,whichCell,clu_rf,best_fit,best_params,r_squared, ...
            rf_stats,rf_rsq,flat_rsq,empirical_rsq,bestmodel_rsq,cluster_info)
        dispstat(['Finished ' thisExp ', progress saved'],'keepprev','keepthis','timestamp')

    end
    wc{k} = whichCell;
    crf{k} = clu_rf;
    bf{k} = best_fit;
    bp{k} = best_params;
    rs{k} = r_squared;
    rfs{k} = rf_stats;
    rfrs{k} = rf_rsq;
    null_rs{k} = flat_rsq;
    emprcl_rs{k} = empirical_rsq;
    best_rs{k} = bestmodel_rsq;
    clinf{k} = cluster_info;
    dispstat(['dataset ' num2str(k) '/' num2str(nDsets) ' complete'],'keepprev','keepthis','timestamp')
end

clu_fits.RF =  vertcat(crf{:});
clu_fits.bestFit = vertcat(bf{:});
clu_fits.bestParams = vertcat(bp{:});
clu_fits.score = vertcat(rs{:});
clu_fits.rfCVscore = vertcat(rfrs{:});
clu_fits.nullScore = vertcat(null_rs{:});
clu_fits.empiricalScore = vertcat(emprcl_rs{:});
clu_fits.bestModelScore = vertcat(best_rs{:});
clu_fits.rfStats = cat(1,rfs{:});
clu_fits.clusterInfo = vertcat(clinf{:});
clu_fits.whichCell = vertcat(wc{:});
clu_fits.whichModel = modelStrings;
clu_fits.whichFunc = modelFuncs;
clu_fits.scoreMetric = params.method;
clu_fits.datasets = db;

save([localRoot '\SC_popstructs\' saveName],'clu_fits')

%% border-reading helper
function [scTop, scBottom] = bord(bfname,snname,spks)
% [scTop, scBottom] = bord(bfname,snname,spks)
%

if exist(bfname,'file')
    brdr = readtable(bfname ,'Delimiter','\t', 'FileType', 'text');
    sc = contains(brdr.acronym,'SCs');
    scupper = brdr.upperBorder(sc);
    sclower = brdr.lowerBorder(sc);
    scTop = max(scupper);
    scBottom = min(sclower);
elseif exist(snname,'file')
    [ctxChan, pagchan] = wheresCortex('foo','snrf',snname);
    
    scTop = spks.ycoords(ctxChan);
    scBottom = spks.ycoords(pagchan);
else
    scTop = [];
    scBottom = [];
end

end

%% saving helper
function parsave(fname,whichCell,clu_rf,best_fit,best_params,r_squared,rf_stats,rf_rsq, ...
    flat_rsq,empirical_rsq,bestmodel_rsq,cluster_info)
% because you can't use 'save' in a parfor loop

save(fname,'whichCell','clu_rf','best_fit','best_params','r_squared', ...
    'rf_stats','rf_rsq','flat_rsq','empirical_rsq','bestmodel_rsq','cluster_info')
end

