% Script that runs through cells to fit multiple RF models and plot the
% results. It is set up to show 3 models, the residuals from the 'best'
% model, and the waveform/stability of the unit.
%
% Requires: - 'snrf' structure made by the getRFs script
%           - 'RFs' folder, which has the models and fitting function

%% params
localRoot = 'C:\DATA\Spikes\';

overwrite = true;
verbose = true;

whichModels = {'flat','gauss','DoG','\nabla gauss','\nabla^2 gauss'}; % model labels
clims = 0.98; %th quantile (+/-), proportion of RF within the colorbar limits

init_std = 1.5; % units of 'stixels', for initial conditions of all models
min_btwn = 2; % stixels; if max and min are closer than this, DoG means start off equal

clear db
ephys_bilateral_db

%%
for k = 1:length(db)
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        thisExp =  [db(k).mouse_name '_' db(k).date '_' thisTag];
        if verbose, disp(['Now: ' thisExp]); end
        
        % load everything for that probe
        [dsetFolders,dataDir] = expDirs(db(k).mouse_name,db(k).date,thisTag,db(k).dataServer);
        ksDir = [dataDir '\sorting\'];
        
        snrfFolder = [localRoot dsetFolders '\'];
        snname = [snrfFolder 'sparse_noise_RFs.mat'];
        if exist(snname,'file')
            snrf = loadVar(snname,'snrf');
        else
            disp([thisExp '      doesn''t have an RF file, nevermind'])
            continue
        end
        
        spks = loadNeuralData(ksDir);
        noiseStart = min(snrf.stimTimes_local);
        noiseEnd = max(snrf.stimTimes_local);
        
        yPos = snrf.YPos; nX = length(yPos); 
        xPos = snrf.XPos; nY = length(xPos);
        y0 = {0,spks.ycoords(max(snrf.responsive_channels))};
        whichY0 = ~isempty(y0{2})+1; % if we can infer sSC top, center spike depth to that
        depthStr = {'(depth from bottom)','(depth from SC surface)'}; % label for posterity
        
        if verbose, disp(' Fitting for all cells ... '); tic; end
        
        fig = figure('units','normalized','position',[0.0901 0.0454 0.7234 0.8685],'visible','off');
        for iClu = 1:length(snrf.neur_ID)
            CID = snrf.neur_ID(iClu);
            
            fname = sprintf('%s_Cell_%d_RFfits.png',thisTag,CID);
            figDir = [localRoot '\' db(k).mouse_name '\' db(k).date '\receptiveFields\'];
            if ~exist(figDir,'dir')
                mkdir(figDir)
            elseif exist([figDir fname],'file') && (~overwrite)
                continue
            end
            
            % fit things
            rf = snrf.neur_rfmap{iClu};
            
            maxZ = (max(rf(:))-mean(rf(:)))./std(rf(:));
            minZ = (min(rf(:))-mean(rf(:)))./std(rf(:));
            flipped = abs(minZ)>maxZ;
            if abs(minZ)>maxZ
                % peak is negative - neuron is suppressed
                useRF = -rf;
            else
                useRF = rf;
            end
            useRF = useRF-quantile(useRF(:),0.25);
            
            [maxY, maxX] = find(useRF==max(useRF(:)),1);
            [minY, minX] = find(useRF==min(useRF(:)),1);
            
            % gaussian fit
            lb = [0,min(xPos),0,min(yPos),0,-pi/4];
            ub = [realmax('double'),max(xPos),(max(xPos))/2,max(yPos),(max(yPos))/2,pi/4];
            
            x0 = [1,xPos(maxX),mean(diff(xPos))*init_std,yPos(maxY),mean(diff(yPos))*init_std,0];
            
            [x, gaus_fit] = fit2dGaussRF_mine(xPos, yPos, useRF,@D2GaussFunctionRot,[x0;lb;ub]);
            gaus_fit = ((-1)^flipped)*gaus_fit;
            
            % difference-of-gaussians
            lb = [0,min(xPos),0,min(yPos),0,0,min(xPos),0,min(yPos),0,-pi/4];
            ub = [realmax('double'),max(xPos),(max(xPos))/2,max(yPos),(max(yPos))/2, ...
                realmax('double'),max(xPos),(max(xPos))/2,max(yPos),(max(yPos))/2,pi/4];
            
            if norm([maxY-minY, maxX-minX]) <= min_btwn
                mu_x1 = mean(xPos([maxX minX])); mu_x2 = mu_x1;
                mu_y1 = mean(yPos([maxY minY])); mu_y2 = mu_y1;
            else
                mu_x1 = xPos(maxX); mu_x2 = xPos(minX);
                mu_y1 = yPos(maxY); mu_y2 = yPos(minY);
            end
            x0 = [1, mu_x1, mean(diff(xPos))*init_std, mu_y1, mean(diff(yPos))*init_std, ...
                1, mu_x2, mean(diff(xPos))*init_std, mu_y2, mean(diff(yPos))*init_std, 0];
            
            [x_dg,dg_fit] = fit2dGaussRF_mine(xPos, yPos, useRF, @D2_DoG_unequal_Rot,[x0;lb;ub]);
            dg_fit = ((-1)^flipped)*dg_fit;
            
            % gradient of a gaussian
            lb = [0, min(xPos), 0, min(yPos), 0, -pi]; % they get more flexibility for rotation
            ub = [realmax('double'), max(xPos), (max(xPos))/2, max(yPos), (max(yPos))/2, pi];
            
            x0 = [1,xPos(maxX),mean(diff(xPos))*init_std,yPos(maxY),mean(diff(yPos))*init_std,pi/2];

            [x_gg,gg_fit] = fit2dGaussRF_mine(xPos, yPos, useRF, @D2_gradGaus_Rot,[x0;lb;ub]);
            gg_fit = ((-1)^flipped)*gg_fit;
            
            % laplacian of a gaussian
            lb = [0, min(xPos), min(diff(xPos))/2, min(yPos), min(diff(yPos))/2, -pi/2];
            ub = [realmax('double'), max(xPos), (max(xPos))/2, max(yPos), (max(yPos))/2, pi/2];
            
            x0 = [1,xPos(maxX),mean(diff(xPos))*init_std,yPos(maxY),mean(diff(yPos))*init_std,0];
            
            [x_lg,lg_fit] = fit2dGaussRF_mine(xPos, yPos, -useRF, @D2_LoG_Rot,[x0;lb;ub]);
            lg_fit = ((-1)^flipped)*lg_fit;
            
            % compare to find 'best' model (for now least-squares, no discounting)
            allmods = {mean(rf(:)), gaus_fit, dg_fit, gg_fit, lg_fit};
            
            H0 = norm(rf(:)-mean(rf(:)))^2; % null model (flat distribution)
            errs = [H0, norm(rf(:)-gaus_fit(:))^2, norm(rf(:)-dg_fit(:))^2, ...
                 norm(rf(:)-gg_fit(:))^2,  norm(rf(:)-lg_fit(:))^2];
            [~,best_model_for_now] = min(errs);
            
            resid = rf - allmods{best_model_for_now}; 
            
            %% make plots 
            hold off;
            
            subtightplot(4,3,1,[0.05,0.05]) % receptive field
            imagesc(xPos,yPos,rf)
            xticks([]); set(gca,'ydir','normal')
            title('receptive field (to white square onset)')
%             hold on; plot([0 0],ylim,'--','color','k')
%             plot(xlim,[0 0],'--','color','k'); hold off
            q = abs(quantile(abs(rf(:)),clims));
            caxis([-q q])
            
            subtightplot(4,3,2,[0.05,0.05]) % timecourse
            plot(snrf.neur_rfstats(iClu).timeBins, snrf.neur_rfstats(iClu).timeCourse, ...
                'color','k','linewidth',2)
            xticks([0 max(xlim)]); box off
            title('time-course from SVD')

            subtightplot(4,3,4,[0.05,0.05]) % model 2 (gaussian)
            imagesc(xPos,yPos,gaus_fit)
            xticks([]); set(gca,'ydir','normal')
            title(['gaussian fit (||err||^2: ' num2str(errs(2),2) ')'])
%             hold on; plot([0 0],ylim,'--','color','k')
%             plot(xlim,[0 0],'--','color','k'); hold off
            q = abs(quantile(abs(gaus_fit(:)),clims));
            caxis([-q q])
            
            subtightplot(4,3,5,[0.05,0.05]) % model 4
            imagesc(xPos,yPos,gg_fit)
            xticks([]); set(gca,'ydir','normal')
            title(['gradient of gauss. (||err||^2: ' num2str(errs(4),2) ')'])
%             hold on; plot([0 0],ylim,'--','color','k')
%             plot(xlim,[0 0],'--','color','k'); hold off
            q = abs(quantile(abs(gg_fit(:)),clims));
            caxis([-q q])
            
            subtightplot(4,3,7,[0.05,0.05]) % model 3
            imagesc(xPos,yPos,dg_fit)
            xticks([]); set(gca,'ydir','normal')
            title(['difference of two gauss. (||err||^2: ' num2str(errs(3),2) ')'])
%             hold on; plot([0 0],ylim,'--','color','k')
%             plot(xlim,[0 0],'--','color','k'); hold off
            q = abs(quantile(abs(dg_fit(:)),clims));
            caxis([-q q])
            
            subtightplot(4,3,8,[0.05,0.05]) % model 5 
            imagesc(xPos,yPos,lg_fit)
            xticks([]); set(gca,'ydir','normal')
            title(['laplacian of gauss. (||err||^2: ' num2str(errs(5),2) ')'])
%             hold on; plot([0 0],ylim,'--','color','k')
%             plot(xlim,[0 0],'--','color','k'); hold off
            q = abs(quantile(abs(lg_fit(:)),clims));
            caxis([-q q])
            
            subtightplot(4,3,10,[0.05,0.05]) % residuals plot
            imagesc(xPos,yPos,resid); set(gca,'ydir','normal')
            title(['residuals from least ^2s model (' whichModels{best_model_for_now} ')'])
%             hold on; plot([0 0],ylim,'--','color','k')
%             plot(xlim,[0 0],'--','color','k'); hold off
            q = abs(quantile(abs(resid(:)),0.98));
            caxis([-q q])
            
            subtightplot(4,3,[11 12],[0.05,0.05]) % amplitude plot
            deezspks = (spks.clu == CID) & (spks.st>noiseStart & spks.st<noiseEnd);
            scatter(spks.st(deezspks),spks.gain*spks.spikeAmps(deezspks),24,'markeredgealpha',0.4)
            set(gca,'yaxislocation','right')
            text(max(spks.st(deezspks)),max(ylim),'spikes during sparse noise', 'fontsize',13, ...
                'horizontalalignment','right','verticalalignment','top')
            
            subtightplot(4,3,[3 6 9],[0.05,0.05]) % waveform
            clut = spks.cluTemps(CID+1)+1;
            [~,ind] = max(max(spks.tempsUnW(clut,:,:)));
            plotchans = rectify(ind-7,1):saturate(ind+7,374);
            wf = permute(spks.tempsUnW(clut,:,plotchans),[3 2 1]);
            plotWaveform(wf,spks.xcoords(plotchans),spks.ycoords(plotchans)-y0{whichY0},20,5,0,'k')
            set(gca,'yaxislocation','right')
            title(sprintf('Cell %d, %s %s %s',CID,db(k).mouse_name,db(k).date,thisTag))
            text(mean(xlim),max(ylim)-5,['Spike template ' depthStr{whichY0}],'fontsize',13, ...
                'horizontalalignment','center','verticalalignment','top')
            
            %% save
            print(fig,[figDir fname],'-dpng')
        end
        close(fig)
        if verbose, disp(['     ... done in ' num2str(toc)]); end
    end
end
if verbose, disp('all done.'); end

%% for later
%             % re-fit on residuals (once I can automatically choose the best model)
%             resid = rf - fit_dg;
%             maxZ = (max(resid(:))-mean(resid(:)))./std(resid(:));
%             minZ = (min(resid(:))-mean(resid(:)))./std(resid(:));
%             if abs(minZ)>maxZ
%                 % peak is negative - neuron is suppressed
%                 useRF_refit = -resid;
%             else
%                 useRF_refit = resid;
%             end
%             useRF_refit = useRF_refit-quantile(useRF_refit(:),0.25);
%
%             [x_res,fit_res] = fit2dGaussRF_mine(yPos, xPos, useRF_refit);
%
%             lb = [0,min(yPos),0,min(xPos),0,0,min(yPos),0,min(xPos),0,-pi/4];
%             ub = [realmax('double'),max(yPos),(max(yPos))^2,max(xPos),(max(xPos))^2, ...
%                 realmax('double'),max(yPos),(max(yPos))^2,max(xPos),(max(xPos))^2,pi/4];
%
%             [maxY, maxX] = find(useRF_refit==max(useRF_refit(:)),1);
%             [minY, minX] = find(useRF_refit==min(useRF_refit(:)),1);
%             x0 = [1,yPos(maxX),mean(diff(yPos))*3,xPos(maxY),mean(diff(xPos))*3, ...
%                 1,yPos(minX),mean(diff(yPos))*3,xPos(minY),mean(diff(xPos))*3,0];
%
%             [x_dg_res,fit_dg_res] = fit2dGaussRF_mine(yPos, xPos, useRF_refit, @D2_DoG_unequal_Rot,[x0;lb;ub]);


%% for even later

% % some other models to try
% lb = [0,min(yPos),0,min(xPos),0,0,min(yPos),min(xPos),-pi/4];
% ub = [realmax('double'),max(yPos),(max(yPos))^2,max(xPos),(max(xPos))^2, ...
%     realmax('double'),max(yPos),max(xPos),pi/4];
% 
% [maxY, maxX] = find(useRF==max(useRF(:)),1);
% [minY, minX] = find(useRF==min(useRF(:)),1);
% x0 = [1,yPos(maxX),mean(diff(yPos))*init_std,xPos(maxY),mean(diff(xPos))*init_std, ...
%     1,minX,minY,0];
% 
% [x_dog,dog_fit] = fit2dGaussRF_mine(yPos, xPos, useRF, @D2_DoG_Rot,[x0;lb;ub]);

% lb = [0,min(yPos),0,min(xPos),0,-pi/4];
% ub = [realmax('double'),max(yPos),(max(yPos))^2,max(xPos),(max(xPos))^2,pi/4];
% 
% [maxY, maxX] = find(useRF==max(useRF(:)),1);
% x0 = [1,yPos(maxX),mean(diff(yPos))*init_std,xPos(maxY),mean(diff(xPos))*init_std,0];
% 
% [x_gg,gg_fit] = fit2dGaussRF_mine(yPos, xPos, useRF, @D2_gradGaus_Rot,[x0;lb;ub]);

% lb = [0,min(yPos),0,min(xPos),0,-pi/4];
% ub = [realmax('double'),max(yPos),(max(yPos))^2,max(xPos),(max(xPos))^2,pi/4];
% 
% [maxY, maxX] = find(useRF==max(useRF(:)),1);
% x0 = [1,yPos(maxX),mean(diff(yPos))*2,xPos(maxY),mean(diff(xPos))*2,0];
% 
% [x_lg,fit_log] = fit2dGaussRF_mine(yPos, xPos, useRF, @D2_LoG_Rot,[x0;lb;ub]);

