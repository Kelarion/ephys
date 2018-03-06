%% Parameters
iSet = 4;
plotDir = 'E:\matteodata\poster\behavior\';
colparts = [1 .2 0.3; % RFL
            1 .6 0; % LFL
            0 .7 1;  % LHL
            .3 0 .7; % RHL
            .1 .9 .6]; % BT

%% 3D trajectories  
close all
niceLims = [123.5733  281.7557];
RFLFig = behTraj(traj_all{iSet}.traj_RFL_smooth, colparts(1,:), niceLims);
LFLFig = behTraj(traj_all{iSet}.traj_LFL_smooth, colparts(2,:), niceLims);
RHLFig = behTraj(traj_all{iSet}.traj_RHL_smooth, colparts(3,:), niceLims);
LHLFig = behTraj(traj_all{iSet}.traj_LHL_smooth, colparts(4,:), niceLims);
BTFig = behTraj(traj_all{iSet}.traj_BT_smooth, colparts(5,:), niceLims);

hgexport(RFLFig,[plotDir '3d_traj_' 'RFL'],...
                hgexport('factorystyle'), 'Format', 'epsc');
hgexport(LFLFig,[plotDir '3d_traj_' 'LFL'],...
                hgexport('factorystyle'), 'Format', 'epsc');
hgexport(RHLFig,[plotDir '3d_traj_' 'RHL'],...
                hgexport('factorystyle'), 'Format', 'epsc');
hgexport(LHLFig,[plotDir '3d_traj_' 'LHL'],...
                hgexport('factorystyle'), 'Format', 'epsc');
hgexport(BTFig,[plotDir '3d_traj_' 'BT'],...
                hgexport('factorystyle'), 'Format', 'epsc');        
            
            hgexport(fig,[plotDir 'nonmossy_squiggles'],...
                hgexport('factorystyle'), 'Format', 'svg');    
            
%% Video
chosenFrame = 16397;
vidFig = figure('Position',[379 532 1109 444]);
for iCam = 0:1
    whichCam = sprintf('bias_video_cam_%d_',iCam);
    trkDir = dir(['\\dm11\hantmanlab\matteo\data\cerebellum_probes\muaddib\6-19-17_2140\' whichCam '*.trk']); % get the track file
    vidDir = dir(['\\dm11\hantmanlab\matteo\data\cerebellum_probes\muaddib\6-19-17_2140\' whichCam '*.avi']); % get the video file
    vidPlace = [vidDir.folder '\' vidDir.name];
    trkPlace = [trkDir.folder '\' trkDir.name];
    
    ptrk = load(trkPlace,'-mat');
    vidObj = VideoReader(vidPlace);

    vidObj.CurrentTime = chosenFrame/vidObj.FrameRate;
    frem = readFrame(vidObj);
    vidAx = subaxis(1,2,(iCam)+1, 'Spacing',0.03, 'Padding', 0, 'Margin', 0.03);
    image(frem)
    hold on;
    for i = 1:size(ptrk.pTrk,1)
        plot(ptrk.pTrk(i,1,chosenFrame),ptrk.pTrk(i,2,chosenFrame),'o', ...
            'Color',colparts(i,:),'MarkerSize', 16*(iCam + 1),'LineWidth',4*(iCam + 1))
    end
    axis(vidAx,'off')
end

figPlace = fullfile(plotDir,[ 'tracked_frame_', num2str(chosenFrame)]);
print(vidFig,figPlace,'-dtiff');

