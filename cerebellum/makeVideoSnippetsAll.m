%% load
load('C:\Users\Matteo A\Documents\MATLAB\mine\ephys\cerebellum\allrips.mat')
ksRoot = 'I:\data\';
localKsRoot = 'C:\DATA\Spikes\Janelia\';
saveDir = 'C:\DATA\Spikes\Janelia\Videos\';

% To convert from the 3-number identification (which correspond to cell ID,
% tag, and dataset) to a single unique cell ID, I use the below mappings. 
% They're similar to a conversion between base 5 and 10. For some 
% reason, using base-5 makes this mapping always invertible. (I checked 
% this with randomly-generated ID matrices where the 'digits' could vary 
% between 1 and 2000.)
b = 5;
wc2id = @(WC) WC*[b^2 b^1 b^0]';
id2wc = @(id) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)];

CID = wc2id(R.whichCell);
cgs = unique(CID);
nClu = length(cgs);

%% assess significance of phase-locking 
% in reality, testing for non-uniformity of angle distribution
minNSpk = 50;
alph = 0.05;

pvals = nan(nClu,3); 
for iClu = 1:nClu
    if (sum(CID == cgs(iClu)) > minNSpk)
        pvals(iClu,1) = circ_rtest(R.spkPhase(CID == cgs(iClu)));
    end
    if (sum(CID == cgs(iClu) & ~R.isOpto) > minNSpk)
        pvals(iClu,2) = circ_rtest(R.spkPhase(CID == cgs(iClu) & ~R.isOpto));
    end
    if (sum(CID == cgs(iClu) & R.isOpto) > minNSpk)
        pvals(iClu,3) = circ_rtest(R.spkPhase(CID == cgs(iClu) & R.isOpto));
    end
end

islocked = sum(pvals(:,[1 2]) < alph,2) == 2; % we care about phase-locking spontaneously
islocked = islocked & (pvals(:,3) < alph | isnan(pvals(:,3)));

%% compile video snippets around spontaneous spiking events
db = R.datasets;
win = [-0.5 0.5];
minBetween = 0; % seconds
overwrite = false;

foo = id2wc(cgs(islocked)); % which datasets contain events
dsets = unique(foo(:,[2 3]),'rows');

for iSet = 1:size(dsets,1)
    k = dsets(iSet,2);
    d = dsets(iSet,1);
    
    vidname = [db(k).name '_' db(k).date '_' db(k).depth{d} '_events.avi'];
    if exist([saveDir '\' vidname],'file') && ~overwrite
        continue 
    end
    
    deez = (R.whichCell(:,2) == d) & (R.whichCell(:,3) == k);
    deez = deez & ~R.isOpto; % must be spontaneous
    deez = deez & ismember(CID,cgs(islocked));
    
    parentDir = [ksRoot '\' db(k).name '\'];
    subDir = dir([parentDir '*' db(k).date]);
    
    if length(subDir) > 1
        subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
        if isempty(subDir)
            disp('can''t find it ¯\_(?)_/¯')
        end
    end
    vidDir = [parentDir subDir.name '\'];
    dataDir = [localKsRoot '\' db(k).name '\' subDir.name '\'];
    
    sp = loadJRCdir(dataDir,true);
    
    vidName = dir([vidDir sprintf('bias_video_cam_%d_',0) '*.avi']);
    vidPlace = [vidName.folder '\' vidName.name]; % get the video files
    vid0 = VideoReader(vidPlace);
    
    vidName = dir([vidDir sprintf('bias_video_cam_%d_',1) '*.avi']);
    vidPlace = [vidName.folder '\' vidName.name]; % get the video files
    vid1 = VideoReader(vidPlace);
    
    if ~isfield(sp,'frameTimes')
        nFrame = round(vid0.Duration*vid0.FrameRate);
        
        % From Britton's code:
        % If there are more frames than pulses, assume that the neural acquisition
        % was terminated before the pulses were stopped, and add timestamps using
        % the average inter-pulse interval.
        i_pulse = sp.i_pulse;
        fs_vid = round(vid0.FrameRate);
        if length(i_pulse) < nFrame
            disp([num2str(nFrame-length(i_pulse)) ' more frames than pulses were counted.']) 
            disp('Assuming neural data was stopped before the movie; extending pulse train ... ')
            i_pulse = [i_pulse i_pulse(end) + median(diff(i_pulse))*(1:(nFrame-length(i_pulse)))];
        end
        % Assume that only initial frames are dropped? DEFINITELY CHECK THIS!!!!!
        t_frame = i_pulse((length(i_pulse)-nFrame+1):end)*1000/sp.sample_rate;
        save([dataDir '\clustered\t_frame.mat'],'t_frame')
        sp.frameTimes = t_frame;
    end
    
    ft = sp.frameTimes/1000;
    vidOffset = ft(1);
    
    % events are detected on single-cell basis, so the same event might 
    % show up on multiple channels. here I merge these times into one 
    % event time. to do this I detect overlapping events as those which
    % start before another has ended, then merge those events, and
    % repeat this until no more events overlap.
    et = unique(R.interval(deez,:),'rows'); % this is sorted by first column (i.e. onset time)
    overlap = find(et(2:end,1) - et(1:end-1,2) < minBetween) + 1;
    while ~isempty(overlap) % we'll take the average across cells
        et(overlap-1,:) = 0.5*(et(overlap-1,:) + et(overlap,:));
        et(overlap,:) = [];  
        overlap = find(et(2:end,1) - et(1:end-1,2) < minBetween) + 1;
    end
    
    %% read and write
    disp(['Now writing ' vidname])
    v = VideoWriter([saveDir  '\' vidname]);
    v.FrameRate = vid0.FrameRate;
    open(v)
    
    figure('units','normalized','position',[0.1240 0.2889 0.6995 0.5250],'visible','off');
    cam0 = subtightplot(1,2,1,[0 0.005],0,0);
    cam1 = subtightplot(1,2,2,[0 0.005],0,0);
    tic
    for iEv = 1:size(et,1)
        t_start = et(iEv,1)-vidOffset;
        t_fin = et(iEv,2)-vidOffset;
        if (t_start<win(1)) || (t_fin > vid0.Duration)
            continue
        end
        vid0.CurrentTime = rectify(t_start+win(1));
        vid1.CurrentTime = rectify(t_start+win(1));
        
        ff = 1;
        while vid0.CurrentTime < min([t_fin+win(2) vid0.Duration])
            vid0Frame = readFrame(vid0);
            vid1Frame = readFrame(vid1);
            if ff == 1 % faster to make the objects only for the first frame
                h0 = image(cam0,vid0Frame);
                h1 = image(cam1,vid1Frame);
                hold(cam1,'on')
                hold(cam0,'on')
                tstamp = text(cam0,vid0.width,vid0.height,num2str(vid0.CurrentTime),...
                    'color','yellow','horizontalalignment','right', ...
                    'verticalalignment','bottom','fontsize',15);
                estamp = text(cam0,1,1,['Event #: ' num2str(iEv)],'color','yellow',...
                    'horizontalalignment','left','verticalalignment','top',...
                    'fontsize',13);
                dot0 = scatter(cam0,400,50,800,'y','filled'); dot0.Visible = 'off';
                dot1 = scatter(cam1,400,50,800,'y','filled'); dot1.Visible = 'off';
                set(cam0,'Visible','off');
                set(cam1,'Visible','off');
                hold(cam0,'off')
                hold(cam1,'off')
            else % and change the data on subsequent frames
                set(h0,'cdata',vid0Frame);
                set(h1,'cdata',vid1Frame);
                tstamp.String = num2str(vid0.CurrentTime);
            end
            
            if (vid0.CurrentTime > t_start) && (vid0.CurrentTime < t_fin)
                dot0.Visible = 'on';
                dot1.Visible = 'on';
            else
                dot0.Visible = 'off';
                dot1.Visible = 'off';
            end
            
            frem = getframe(gcf);
            writeVideo(v,frem);
            ff = ff+1;
            %     pause(1/vidObj.FrameRate);
        end
        if ~mod(iEv,5) || (iEv == size(et,1))
            disp(['event ' num2str(iEv) ' written'])
        end
    end
    close(gcf)
    close(v)
    disp(['Finished in ' num2str(toc) ' seconds'])
end



