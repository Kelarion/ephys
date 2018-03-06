% Test of full-body tracking and basal pons cell analysis. Objects / cameras
% are:
%
% Cam0: bottom-up view of mouse.
%   [1] Right forelimb
%   [2] Left forelimb
%   [3] Left hindlimb
%   [4] Right hindlimb
%   [5] Base of tail
%
% Cam1: front-on view of mouse.
%   [1] Right forelimb
%   [2] Left forelimb
%
% Cam2: back view of mouse.
%   [1] Left hindlimb
%   [2] Right hindlimb
%   [3] Base of tail
%
% Output dimensions are in the coordinates of cam1:
% [1] 
% [2]
% [3]

%% load data

t_spk = t_spk_all{dset2plot};

t_frame = t_frame_all{dset2plot};

traj_RFL = traj_all{dset2plot}.traj_RFL;
traj_LFL = traj_all{dset2plot}.traj_LFL;
traj_RHL = traj_all{dset2plot}.traj_RHL;
traj_LHL = traj_all{dset2plot}.traj_LHL;
traj_BT = traj_all{dset2plot}.traj_BT;
traj_RFL_smooth = traj_all{dset2plot}.traj_RFL_smooth;
traj_LFL_smooth = traj_all{dset2plot}.traj_LFL_smooth;
traj_RHL_smooth = traj_all{dset2plot}.traj_RHL_smooth;
traj_LHL_smooth = traj_all{dset2plot}.traj_LHL_smooth;
traj_BT_smooth = traj_all{dset2plot}.traj_BT_smooth;

fs_vid = round(1000/mean(diff(t_frame_all{dset2plot}([1 find(diff(t_frame_all{dset2plot}) < 100)]))));

% Speed and velocity.
v_RFL = [[0 0 0]' diff(traj_RFL_smooth,1,2)]*fs_vid;
v_LFL = [[0 0 0]' diff(traj_LFL_smooth,1,2)]*fs_vid;
v_RHL = [[0 0 0]' diff(traj_RHL_smooth,1,2)]*fs_vid;
v_LHL = [[0 0 0]' diff(traj_LHL_smooth,1,2)]*fs_vid;
v_BT = [[0 0 0]' diff(traj_BT_smooth,1,2)]*fs_vid;
spd_RFL = sqrt(sum(v_RFL.^2,1));
spd_LFL = sqrt(sum(v_LFL.^2,1));
spd_RHL = sqrt(sum(v_RHL.^2,1));
spd_LHL = sqrt(sum(v_LHL.^2,1));
spd_BT = sqrt(sum(v_BT.^2,1));

%% Plot position vs time.
space_pos = 25;
raw_plot_lims = [4 5]*60*30000; % Can't plot all the raw data at once, so just choose a segment.
for iii = 1:length(traj_all)
    figure('position',[0 0 1024 768],'Visible','off')
    
    h_timeseries{iii}(1) = subplot(6,1,1);
    for j = 1:3;
        plot(t_frame_all{iii},traj_all{iii}.traj_RFL_smooth(j,:)-mean(traj_all{iii}.traj_RFL_smooth(j,:))...
            -(j-1)*space_pos,'color',col_parts(1,:)); hold on
    end
    
    h_timeseries{iii}(2) = subplot(6,1,2);
    for j = 1:3;
        plot(t_frame_all{iii},traj_all{iii}.traj_LFL_smooth(j,:)-mean(traj_all{iii}.traj_LFL_smooth(j,:))...
            -(j-1)*space_pos,'color',col_parts(2,:)); hold on
    end
    
    h_timeseries{iii}(3) = subplot(6,1,3);
    for j = 1:3;
        plot(t_frame_all{iii},traj_all{iii}.traj_RHL_smooth(j,:)-mean(traj_all{iii}.traj_RHL_smooth(j,:))...
            -(j-1)*space_pos,'color',col_parts(3,:)); hold on
    end
    
    h_timeseries{iii}(4) = subplot(6,1,4);
    for j = 1:3;
        plot(t_frame_all{iii},traj_all{iii}.traj_LHL_smooth(j,:)-mean(traj_all{iii}.traj_LHL_smooth(j,:))...
            -(j-1)*space_pos,'color',col_parts(4,:)); hold on
    end
    
    h_timeseries{iii}(5) = subplot(6,1,5);
    for j = 1:3;
        plot(t_frame_all{iii},traj_all{iii}.traj_BT_smooth(j,:)-mean(traj_all{iii}.traj_BT_smooth(j,:))...
            -(j-1)*space_pos,'color',col_parts(5,:)); hold on
    end
    
    
    h_timeseries{iii}(6) = subplot(6,1,6);
    %{
    plot_SS_CS(neuro_files_all{iii},raw_plot_lims)
    ylim([-.006 .006])
    %}
    raster(t_spk); hold on
    
end
for iii = 1:length(traj_all)
    linkaxes(h_timeseries{iii},'x')
    for i = 1:5
        set(h_timeseries{iii}(i),'position',[0 1-i*.12 1 .12],'xtick',[],'ytick',[])
        axes(h_timeseries{iii}(i))
        set(gca,'ylim',[-3*space_pos space_pos])
        set(zoom,'Motion','horizontal','Enable','on');
        set(pan,'Motion','horizontal','Enable','on')
    end
    set(h_timeseries{iii}(6),'position',[0 0 1 .4],'xtick',[],'ytick',[])
end


%% Regression of simple spike counts on velocity

for i = 1:length(traj_all)
    nf = neuro_files_all{i};
    slashes = strfind(nf,'\');
    neurodir = nf(1:slashes(end));
    plotdir = [neurodir 'plots\regression\'];
    for j = 1:length(t_spk_all{i})
        close all
        deez_spks = t_spk_all{i}(j);
        if all(deez_spks{1} < t_frame_all{i}) % make sure the cell has spikes during the movie
            fprintf('Oops! looks like cell %d doesn''t spike during the movie.\n', j)
            fprintf('Moving on ... \n')
            continue
        end
        [beta_all{i}{j} beta_int_all{i}{j} stats_out_all{i}{j}] = ...
            matteo_regression(traj_all{i},t_spk_all{i}(j),t_frame_all{i},-1000:10:1000,50);
        
        figure(3)
        fig_hand = gcf;
        fig_hand.Color = [1 1 1];
        fig_hand.InvertHardcopy = 'off';
        set(fig_hand,'PaperPositionMode','auto')
        hgexport(fig_hand,[plotdir 'betas_3D_cell_' num2str(j)],...
                hgexport('factorystyle'), 'Format', 'epsc');
            
        figure(2)
        fig_hand = gcf;
        fig_hand.Color = [1 1 1];
        fig_hand.InvertHardcopy = 'off';
        set(fig_hand,'PaperPositionMode','auto')
        hgexport(fig_hand,[plotdir 'betas_1D_cell_' num2str(j)],...
                hgexport('factorystyle'), 'Format', 'epsc');
            
        figure(1)
        fig_hand = gcf;
        fig_hand.Color = [1 1 1];
        fig_hand.InvertHardcopy = 'off';
        set(fig_hand,'PaperPositionMode','auto')
        hgexport(fig_hand,[plotdir 'R2_cell_' num2str(j)],...
                hgexport('factorystyle'), 'Format', 'epsc');
            
        display(num2str(j));
    end
end

% % Summary of tuning for all neurons and effectors.
% for i = 1:length(beta_all)
%     [max_R2(i) i_max_R2(i)] = max(stats_out_all{i}(:,1));
%     beta_max(i,:) = beta_all{i}(i_max_R2(i),2:end);
% end

% Plot summary.
% figure('position',[0 0 450*[5.5 2]],'Visible','off')
% %mylim = [-.0125 .0125];
% for i = 1:length(beta_all)
%     h_3d(i) = subplot(2,5,i);
%     
%     mylim = max(abs(beta_max(i,:)))*[-1.4 1.4];
%     x_circ = .015*diff(mylim)*cos(0:pi/100:2*pi);
%     y_circ = .015*diff(mylim)*sin(0:pi/100:2*pi);
%     
%     plot3([-2 2],[0 0],mylim(1)*[1 1],'--','color',[.8 .8 .8]); hold on
%     plot3([0 0],[-2 2],mylim(1)*[1 1],'--','color',[.8 .8 .8]); hold on
%     plot3([-2 2],mylim(2)*[1 1],[0 0],'--','color',[.8 .8 .8]); hold on
%     plot3([0 0],mylim(2)*[1 1],[-2 2],'--','color',[.8 .8 .8]); hold on
%     plot3(mylim(2)*[1 1],[-2 2],[0 0],'--','color',[.8 .8 .8]); hold on
%     plot3(mylim(2)*[1 1],[0 0],[-2 2],'--','color',[.8 .8 .8]); hold on
% 
%     for ii = 1:5
%         plot3([0 beta_max(i,(ii-1)*3+1)],[0 beta_max(i,(ii-1)*3+2)],[0 beta_max(i,(ii-1)*3+3)],'color',col_parts(ii,:),'linewidth',2); hold on
%         plot3(beta_max(i,(ii-1)*3+1),beta_max(i,(ii-1)*3+2),beta_max(i,(ii-1)*3+3),'ok','markerfacecolor',col_parts(ii,:),'markersize',3); hold on
%         
%         plot3(mylim(2)*[1 1],[0 beta_max(i,(ii-1)*3+2)],[0 beta_max(i,(ii-1)*3+3)],'-','color',col_parts_light(ii,:)); hold on
%         plot3([0 beta_max(i,(ii-1)*3+1)],mylim(2)*[1 1],[0 beta_max(i,(ii-1)*3+3)],'-','color',col_parts_light(ii,:)); hold on
%         plot3([0 beta_max(i,(ii-1)*3+1)],[0 beta_max(i,(ii-1)*3+2)],mylim(1)*[1 1],'-','color',col_parts_light(ii,:)); hold on
%         
%         fill3(mylim(2)*ones(1,length(x_circ)),x_circ+beta_max(i,(ii-1)*3+2),y_circ+beta_max(i,(ii-1)*3+3),col_parts_light(ii,:),'edgecolor',[.4 .4 .4]); hold on
%         fill3(x_circ+beta_max(i,(ii-1)*3+1),mylim(2)*ones(1,length(x_circ)),y_circ+beta_max(i,(ii-1)*3+3),col_parts_light(ii,:),'edgecolor',[.4 .4 .4]); hold on
%         fill3(x_circ+beta_max(i,(ii-1)*3+1),y_circ+beta_max(i,(ii-1)*3+2),mylim(1)*ones(1,length(x_circ)),col_parts_light(ii,:),'edgecolor',[.4 .4 .4]); hold on
%     end
% 
%     xlim(mylim); ylim(mylim); zlim(mylim)
%     title(label_dset{i})
%     view([-50 24])
%     xlabel('Right (mm/s)'); ylabel('Forward (mm/s)'); zlabel('Up (mm/s)')
% end

%% Attempt to replicate Streng, Popova, and Ebner (2017): regression on segments of data around CS.
% for i = 1:length(traj_all);
%     [beta_peri_CS_all{i} beta_peri_CS_int_all{i} stats_out_peri_CS_all{i}] = ...
%         regression_velocity_streng(traj_all{i},t_SS_all{i},t_CS_all{i},t_frame_all{i},-100:1:100,100,-1000:10:1000); display(i)
% end

%% Compare distribution of CS-centered features to full distribution.
% t_offset_CS = -2000:10:2000;
% for i = 1:length(traj_all)
%     display(i);
%     
%     peri_CS{i}.v_RFL = zeros(3,length(t_CS_all{i}),length(t_offset_CS));
%     peri_CS{i}.v_LFL = zeros(3,length(t_CS_all{i}),length(t_offset_CS));
%     peri_CS{i}.v_RHL = zeros(3,length(t_CS_all{i}),length(t_offset_CS));
%     peri_CS{i}.v_LHL = zeros(3,length(t_CS_all{i}),length(t_offset_CS));
%     peri_CS{i}.v_BT = zeros(3,length(t_CS_all{i}),length(t_offset_CS));
%     
%     for j = 1:length(t_CS_all{i})
%         for k = 1:3;
%             peri_CS{i}.v_RFL(k,j,:) = interp1(t_frame_all{i},traj_all{i}.v_RFL(k,:),t_CS_all{i}(j)+t_offset_CS);
%             peri_CS{i}.v_LFL(k,j,:) = interp1(t_frame_all{i},traj_all{i}.v_LFL(k,:),t_CS_all{i}(j)+t_offset_CS);
%             peri_CS{i}.v_RHL(k,j,:) = interp1(t_frame_all{i},traj_all{i}.v_RHL(k,:),t_CS_all{i}(j)+t_offset_CS);
%             peri_CS{i}.v_LHL(k,j,:) = interp1(t_frame_all{i},traj_all{i}.v_LHL(k,:),t_CS_all{i}(j)+t_offset_CS);
%             peri_CS{i}.v_BT(k,j,:) = interp1(t_frame_all{i},traj_all{i}.v_BT(k,:),t_CS_all{i}(j)+t_offset_CS);
%         end
%     end
%       
% end
% 
% % Plot.
% for i = 1:length(traj_all)
%     figure('position',[0 0 400*[8 3]],'Visible','off')
%     for k = 1:3
%         subplot(3,5,sub2ind([5 3],1,k))
%         plot(t_offset_CS,nanmean(squeeze(peri_CS{i}.v_RFL(k,:,:)),1),'color',col_parts(1,:));
%         
%         subplot(3,5,sub2ind([5 3],2,k))
%         plot(t_offset_CS,nanmean(squeeze(peri_CS{i}.v_LFL(k,:,:)),1),'color',col_parts(2,:));
%         
%         subplot(3,5,sub2ind([5 3],3,k))
%         plot(t_offset_CS,nanmean(squeeze(peri_CS{i}.v_RHL(k,:,:)),1),'color',col_parts(3,:));
%         
%         subplot(3,5,sub2ind([5 3],4,k))
%         plot(t_offset_CS,nanmean(squeeze(peri_CS{i}.v_LHL(k,:,:)),1),'color',col_parts(4,:));
%         
%         subplot(3,5,sub2ind([5 3],5,k))
%         plot(t_offset_CS,nanmean(squeeze(peri_CS{i}.v_BT(k,:,:)),1),'color',col_parts(5,:));
%     end
% end