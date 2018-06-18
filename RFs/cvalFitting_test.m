localRoot = 'C:\DATA\Spikes\';
modelStrings = {'gaussian','gabor','DoG','DoUG'};
modelFuncs = {@D2GaussFunctionRot, @D2GaborRot, @D2_DoG_Rot, @D2_DoG_unequal_Rot};
ipsiFits = {'gaussian'};% 'same'};

params.nFold = 4;
params.nStoch = 20; 
params.method = 'direct';

ephys_bilateral_db

min_for_ipsi = 40;

testclu = [1149 2 4;
%     1154 2 4;
    66   1 4;
    1089 1 4;
    307  2 6;
    994  2 6;
    1072 3 6;
%     1077 3 6;
    11   3 7;
    1171 1 7;
    351  1 7;
    562  1 8;
%     1104 1 8;
%     1109 1 8;
    1111 1 8;
    1133 1 8;];

clu_rsqr = [];
clu_fits = struct;
%%
prevk = NaN;
prevt = NaN;
for ii = 1:size(testclu,1)
    k = testclu(ii,3);
    t = testclu(ii,2);
    
    if (t ~= prevt) || (k~=prevk)

        thisTag = db(k).tags{t};
        
        thisExp =  [db(k).mouse_name '_' db(k).date '_' thisTag];
%         if verbose, disp(['Now: ' thisExp]); end
        
        % load everything for that probe
        [dsetFolders,dataDir] = expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [dataDir '\sorting\'];
        
        snrfFolder = [localRoot dsetFolders '\'];
        snname = [snrfFolder 'sparse_noise_RFs.mat'];
        if exist(snname,'file')
            snrf = loadVar(snname,'snrf');
        else
            disp([thisExp '      doesn''t have an RF file, nevermind'])
            %     continue
        end
        
        spks = loadNeuralData(ksDir);
        noiseStart = min(snrf.stimTimes_local);
        noiseEnd = max(snrf.stimTimes_local);
        
        yPos = snrf.YPos; nX = length(yPos);
        xPos = snrf.XPos; nY = length(xPos);
    end
    
    prevk = k;
    prevt = t;
     %%
    CID = testclu(ii,1);
    
    iClu = find(snrf.neur_ID == CID);
    
    rf = snrf.neur_rfmap{iClu};
    
    deezspks = (spks.clu == CID) & (spks.st>noiseStart & spks.st<noiseEnd);
    
    nMod = length(modelStrings);
    rsq = nan(1,nMod,params.nFold);
    mdl = struct('String',{},'func',{},'fitIpsi',{},'ipsitype',{});
    mnX0 = zeros(1,nMod);
    mnWx = zeros(1,nMod);
    dispstat('','init')
    dispstat(['Starting cell: ' num2str(CID)],'keepthis','timestamp')
    for iMod = 1:nMod
        mdl(iMod).String = modelStrings{iMod};
        mdl(iMod).func = modelFuncs{iMod};
        mdl(iMod).fitIpsi = false;
        mdl(iMod).ipsitype = 'no';
        
        [rsq(1,iMod,:),x] = cvalGausModel(mdl(iMod),spks.st(deezspks),snrf.stimTimes_local,snrf.stimPosition,params);
        mnX0(iMod) = median(round(x(2,:)));
        mnWx(iMod) = median(round(x(3,:)));
        
        dispstat(['Done with model: ' modelStrings{iMod}],'timestamp')
    end
    
    straw.String = 'rf';
    rf_rsq = mean(cvalGausModel(straw,spks.st(deezspks),snrf.stimTimes_local,snrf.stimPosition,params));
    
    [~,best_model_for_now] = max(mean(rsq,3));
    % don't bother fitting ipsilateral field if the best fit was too medial
    if abs(mnX0(best_model_for_now)) >= min_for_ipsi
        for jj = 1:length(ipsiFits)
            mdl(nMod+jj).String = modelStrings{best_model_for_now};
            mdl(nMod+jj).func = modelFuncs{best_model_for_now};
            mdl(nMod+jj).fitIpsi = true;
            mdl(nMod+jj).ipsitype = ipsiFits{jj};
            
            rsq(1,iMod+jj,:) = cvalGausModel(mdl(nMod+jj),spks.st(deezspks),snrf.stimTimes_local,snrf.stimPosition,params);
            
            dispstat(['Done with model: ' modelStrings{best_model_for_now} ' with ' ipsiFits{jj} ' ipsilaterally'],'timestamp')
            
        end
    else
        rsq(1,iMod+1:iMod+length(ipsiFits),:) = NaN;
    end
    % get the best fits
    [~, absolute_best] = max(mean(rsq,3));
    wantedModels = unique([best_model_for_now, absolute_best]);
    
    [all_fits, prms] = fitTheseModRF(mdl(2),rf,xPos,yPos,1.5,params.nStoch);
%     [all_fits, prms] = fitTheseModRF(mdl(wantedModels),rf,xPos,yPos,1.5,params.nStoch);
    
    clu_fits(ii).all_fits = all_fits;
    clu_fits(ii).r_squared = rsq;
    clu_fits(ii).parameters = prms;
    clu_fits(ii).rf = rf;
    clu_fits(ii).null_rsq = rf_rsq;
    clu_rsqr = cat(1,clu_rsqr,rsq);
    
end
%%
% % figure;
for ii = 1:16
    k = testclu(ii,3);
    t = testclu(ii,2);
    CID = testclu(ii,1);
    
    [~,best_model_for_now] = max(mean(clu_fits(ii).r_squared(1,1:nMod,:),3));
    [~,absolute_best] = max(mean(clu_fits(ii).r_squared,3));
    
    subplot(4,2,1)
    imagesc(xPos,yPos,clu_fits(ii).rf)
    title(sprintf('Cell %d, %s %s %s',CID,db(k).mouse_name,db(k).date,db(k).tags{t}))
    set(gca,'ydir','normal')
    q = abs(quantile(abs(clu_fits(ii).rf(:)),0.95));
    caxis([-q q])
    
    subplot(4,2,2)
    res = clu_fits(ii).all_fits(:,:,absolute_best)-clu_fits(ii).rf;
    imagesc(xPos,yPos,res)
    title(['Residuals from best'])
    set(gca,'ydir','normal')
    q = abs(quantile(abs(res(:)),0.95));
    caxis([-q q])
    
    for mood = 1:nMod
        subplot(4,2,2+mood)
        thisFit = clu_fits(ii).all_fits(:,:,mood);
        imagesc(xPos,yPos,thisFit)
        dr = mean(clu_fits(ii).r_squared(1,mood,:))*100;% - clu_fits(ii).null_rsq ;
        title([modelStrings{mood} ' %var: ' num2str(dr)])
        set(gca,'ydir','normal')
        q = abs(quantile(abs(thisFit(:)),0.95));
        caxis([-q q])
    end
    
    if size(clu_fits(ii).all_fits,3) > nMod
        subplot(4,2,7)
        thisFit = clu_fits(ii).all_fits(:,:,end);
        imagesc(xPos,yPos,thisFit)
        title([modelStrings{best_model_for_now} ' with ipsi gaussian, %var: '  num2str(mean(clu_fits(ii).r_squared(1,end,:))*100)]) % - clu_fits(ii).null_rsq
        set(gca,'ydir','normal')
        q = abs(quantile(abs(thisFit(:)),0.95));
        caxis([-q q])
    else
        subplot(4,2,7)
        imagesc(0)
    end
    
%     subplot(4,2,8)
%     imagesc(xPos,yPos,clu_fits(ii).all_fits(:,:,5))
%     title(['gaus ipsi, R2: ' num2str(clu_fits(ii).r_squared(5)*100 )])
%     set(gca,'ydir','normal')
    
%     subplot(4,2,3)
%     imagesc(xPos,yPos,rglr_mod)
%     title([modelStrings{foo}])
%     
%     subplot(4,2,4)
%     imagesc(xPos,yPos,rglr_mod + ipsi_mod) 
%     title(['With ' mdl(useThis).ipsitype ' ipsi, R^2: ' num2str(clu_rsqr(ii,useThis) - clu_params(ii).rf_score)])
    
    pause
    
end
