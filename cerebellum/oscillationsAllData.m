%% load 
% load('C:\Users\Matteo A\Documents\MATLAB\mine\ephys\cerebellum\allrips.mat')
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

%% mean phase comparisons
[cidOp,~,indOp] = unique(CID(R.isOpto));
[optoPhs,optoMag] = splitapply(@circmean,R.spkPhase(R.isOpto),indOp);

[cidSp,~,indSp] = unique(CID(~R.isOpto));
[spontPhs,spontMag] = splitapply(@circmean,R.spkPhase(~R.isOpto),indSp);
spontFreq = splitapply(@mean,R.freq(~R.isOpto),indSp);
spontPk = splitapply(@mean,R.magnitude(~R.isOpto),indSp);

inboth = intersect(cidSp,cidOp); % which cells fired in both events
spNop = ismember(cidSp,inboth);
opNsp = ismember(cidOp,inboth); % get indices of these cells
spontPhs = spontPhs(spNop);
spontMag = spontMag(spNop);
spontFreq = spontFreq(spNop);
spontPk = spontPk(spNop);
optoPhs = optoPhs(opNsp);
optoMag = optoMag(opNsp);

lckd = ismember(inboth,cgs(islocked));

mags = spontMag; % determine the magnitude of each point

dtheta = spontPhs - optoPhs; % plot 
phsDiff = mags.*exp(1i*dtheta);
phsors = reshape([zeros(size(phsDiff)) phsDiff nan(size(phsDiff))].',length(inboth)*3,1);

% plot
figure
subplot(1,2,1)
scatlab(optoPhs(lckd),spontPhs(lckd),inboth(lckd),64*mags(lckd),spontFreq(lckd),'filled')
hold on
scatlab(optoPhs(~lckd),spontPhs(~lckd),inboth(~lckd),48*mags(~lckd),spontFreq(~lckd))
samelims
xlabel('phase during induced')
ylabel('phase during spontaneous')
title('Mean phase-locking')
grid on; box on
plot(xlim,xlim,'--','color','k')
plot(xlim,xlim - 2*pi,'--','color','k')
plot(xlim,xlim + 2*pi,'--','color','k')
hold off

subplot(1,2,2)
set(gca,'XLim',[-1, 1])
set(gca,'YLim',[-1, 1])
hold on
% make circular grid
for rr = [0.25 0.5 0.75 1]
    plot(rr*sin(0:0.1:2.1*pi),rr*cos(0:0.1:2.1*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
for rr = [0.25 0.5 0.75 1]'
    plot(rr*sin(0:pi/8:2*pi),rr*cos(0:pi/8:2*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
% ------------------
plot(real(phsors),imag(phsors),'linewidth',2)
scatlab(real(phsDiff(lckd)),imag(phsDiff(lckd)),inboth(lckd),[],spontFreq(lckd),'filled')
scatlab(real(phsDiff(~lckd)),imag(phsDiff(~lckd)),inboth(~lckd),[],spontFreq(~lckd))

title('\boldmath${\bar{z}_{spont} / \hat{\bar{z}}_{opto}}$','fontsize',14,'interpreter','latex')
xlabel('Real part');ylabel('Imaginary part')

c = colorbar;
c.Label.String = 'mean spontaneous frequency (Hz)';
c.Label.FontSize = 11;
% colormap jet

%% -----------------------------
foo = cgs(islocked);
eg_cid = foo(randi(length(foo)));

wrapped_spont = exp(1i*(R.spkPhase(CID == eg_cid & ~R.isOpto)));
[theta, r] = circmean(R.spkPhase(CID == eg_cid & ~R.isOpto));
z_spont = r*exp(1i*theta);

wrapped_opto = exp(1i*(R.spkPhase(CID == eg_cid & R.isOpto)));
[theta, r] = circmean(R.spkPhase(CID == eg_cid & R.isOpto));
z_opto = r*exp(1i*theta);

figure
subplot(1,2,1);
xlim([-1 1]);ylim([-1 1]);
for rr = [0.25 0.5 0.75 1]
    plot(rr*sin(0:0.1:2.1*pi),rr*cos(0:0.1:2.1*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
    hold on
end
for rr = [0.25 0.5 0.75 1]'
    plot(rr*sin(0:pi/8:2*pi),rr*cos(0:pi/8:2*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
scatter(real(wrapped_spont),imag(wrapped_spont),'.')
plot([0 real(z_spont)],[0 imag(z_spont)],'linewidth',2)
scatter(real(z_spont),imag(z_spont),'r','filled')
% axis off
hold off

subplot(1,2,2);
xlim([-1 1]);ylim([-1 1]);
for rr = [0.25 0.5 0.75 1]
plot(rr*sin(0:0.1:2.1*pi),rr*cos(0:0.1:2.1*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
hold on
end
for rr = [0.25 0.5 0.75 1]'
plot(rr*sin(0:pi/8:2*pi),rr*cos(0:pi/8:2*pi),'color',[0.8 0.8 0.8],'linewidth',0.1)
end
scatter(real(wrapped_opto),imag(wrapped_opto),'.')
plot([0 real(z_opto)],[0 imag(z_opto)],'linewidth',2)
scatter(real(z_opto),imag(z_opto),'r','filled')
% axis off
hold off
pause


%% timing of events
interval_bins = linspace(0,10,51); % inter-event intervals
oid = unique(wc2id(R.osc.whichCell));

start2start = zeros(length(interval_bins)-1,length(oid));
end2end = zeros(length(interval_bins)-1,length(oid));
end2start = zeros(length(interval_bins)-1,length(oid));

for iEv = 1:length(oid)
    evTimes = R.osc.Times(wc2id(R.osc.whichCell) == oid(iEv) & isnan(R.osc.optoFreq),:);
    start2start(:,iEv) = histcounts(diff(evTimes(:,1)),interval_bins);
    end2end(:,iEv) = histcounts(diff(evTimes(:,2)),interval_bins);
    end2start(:,iEv) = histcounts(diff([evTimes(2:end,1) evTimes(1:end-1,2)]),interval_bins);
end

% % for only looking at the phase-locked events
% start2start = zeros(length(interval_bins)-1,sum(islocked));
% end2end = zeros(length(interval_bins)-1,sum(islocked));
% end2start = zeros(length(interval_bins)-1,sum(islocked));
% 
% cids = cgs(islocked);
% for iClu = 1:sum(islocked)
%     cluTimes = R.interval(CID == cids(iClu) & ~R.isOpto,:);
%     [~,theseEvs] = unique(R.whichEvent(CID == cids(iClu) & ~R.isOpto));
%     start2start(:,iClu) = histcounts(diff(cluTimes(theseEvs,1)),interval_bins);
%     end2end(:,iClu) = histcounts(diff(cluTimes(theseEvs,2)),interval_bins);
%     foo = cluTimes(theseEvs,:);
%     end2start(:,iClu) = histcounts(diff([foo(2:end,1) foo(1:end-1,2)]),interval_bins);
% end

figure;
bar(interval_bins(1:end-1),sum(end2start,2),'histc')
xlabel('Inter-event interval (sec)')
ylabel('Number of events')

% --------------------

durbins = linspace(4,45,41); % spontaneous event duration distribution
n = histcounts(diff(R.osc.Times(isnan(R.osc.optoFreq),:),[],2).*R.osc.Freq(isnan(R.osc.optoFreq)),durbins,...
    'normalization','probability');

durbins = durbins(1:end-1);
linmod = polyfit(durbins(n>0),log10(n(n>0)),1);

figure;
subplot(2,1,1)
bar(durbins,n,'histc')
ylabel('Probability')
xticks([])

subplot(2,1,2)
scatter(durbins,n,'filled')
hold on
plot(xlim,10.^([xlim; 1 1]'*linmod'))
set(gca,'yscale','log')
ylabel('log10(Probability)')
xlabel('Number of cycles in event')
text(max(xlim),max(ylim),['p(cycle) = ' num2str(10^linmod(1))], ... 
    'horizontalalignment','right','verticalalignment','top')

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
            disp('can''t find it �\_(?)_/�')
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

%% Oscillation-triggered average videos
db = R.datasets;
win = [-0.5 0.5]; % time before and after lock point
lockTo = 1; % 1 for onset, 2 for offset
overwrite = false;

foo = id2wc(cgs(islocked)); % which datasets contain events
dsets = unique(foo(:,[2 3]),'rows');

for iSet = 1:size(dsets,1)
    k = dsets(iSet,2);
    d = dsets(iSet,1);
    
    switch lockTo
        case 1
            sufx = '_meanEventOnset';
        case 2
            sufx = '_meanEventOffset';
    end
    vidname = [db(k).name '_' db(k).date '_' db(k).depth{d} sufx '.avi'];
    if exist([saveDir '\' vidname],'file') && ~overwrite
        continue 
    end
    
    deez = (R.whichCell(:,2) == d) & (R.whichCell(:,3) == k);
    deez = deez & ~R.isOpto; % must be spontaneous
    deez = deez & ismember(CID,cgs(islocked)); % and have phase-locked spiking
    deez = find(deez);
    
    parentDir = [ksRoot '\' db(k).name '\'];
    subDir = dir([parentDir '*' db(k).date]);
    
    if length(subDir) > 1
        subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
        if isempty(subDir)
            disp('can''t find it �\_(?)_/�')
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
    
    [et,inds] = unique(R.interval(deez,:),'rows'); % this is sorted by first column (i.e. onset time)
    whichClu = CID(deez(inds));
    clu = unique(whichClu);
    
    %% read and write
    frem_offsets = [win(1):(1/vid0.FrameRate):win(2)];
    
    disp(['Now writing ' vidname])
    v = VideoWriter([saveDir  '\' vidname]);
    v.FrameRate = vid0.FrameRate;
    open(v)
    
    figure('units','normalized','position',[0.1240 0.2889 0.6995 0.5250],'visible','off');
    cam0 = subtightplot(1,2,1,[0 0.005],0,0);
    cam1 = subtightplot(1,2,2,[0 0.005],0,0);
    tic
    for iClu = 1:length(clu)
        cluEvs = et(whichClu == clu(iClu),:);
        for ff = 1:length(frem_offsets)
            avgVid0_frame = zeros(vid0.Height,vid0.Width,3,'uint8');
            avgVid1_frame = zeros(vid1.Height,vid1.Width,3,'uint8');
            for iEv = 1:size(cluEvs,1) % get average over all events
                t_frame = cluEvs(iEv,lockTo) + frem_offsets(ff) - vidOffset;
                if (t_frame<0) || (t_frame > vid0.Duration)
                    continue
                end
                vid0.CurrentTime = t_frame;
                vid1.CurrentTime = t_frame;
                
                vid0Frame = readFrame(vid0);
                vid1Frame = readFrame(vid1);
                
                avgVid0_frame = avgVid0_frame + (1/size(cluEvs,1))*vid0Frame;
                avgVid1_frame = avgVid1_frame + (1/size(cluEvs,1))*vid1Frame;
            end
            
            if ff == 1 % faster to make the objects only for the first frame
                h0 = image(cam0,avgVid0_frame);
                h1 = image(cam1,avgVid1_frame);
                hold(cam1,'on')
                hold(cam0,'on')
                tstamp = text(cam0,vid0.width,vid0.height,num2str(frem_offsets(ff)),...
                    'color','yellow','horizontalalignment','right', ...
                    'verticalalignment','bottom','fontsize',15);
                estamp = text(cam0,1,1,['Neuron #: ' num2str(clu(iClu))],'color','yellow',...
                    'horizontalalignment','left','verticalalignment','top',...
                    'fontsize',13);
                dot0 = scatter(cam0,400,50,800,'y','filled'); dot0.Visible = 'off';
                dot1 = scatter(cam1,400,50,800,'y','filled'); dot1.Visible = 'off';
                set(cam0,'Visible','off');
                set(cam1,'Visible','off');
                hold(cam0,'off')
                hold(cam1,'off')
            else % and change the data on subsequent frames
                set(h0,'cdata',avgVid0_frame);
                set(h1,'cdata',avgVid1_frame);
                tstamp.String = num2str(frem_offsets(ff));
            end
            
            if frem_offsets(ff) >=0 
                dot0.Visible = 'on';
                dot1.Visible = 'on';
            else
                dot0.Visible = 'off';
                dot1.Visible = 'off';
            end
            
            frem = getframe(gcf);
            writeVideo(v,frem);
            %     pause(1/vidObj.FrameRate);
           
        end
        disp(['Cell ' num2str(iClu) ' written at ' num2str(toc) ' seconds'])
    end
    close(gcf)
    close(v)
    disp(['Finished in ' num2str(toc) ' seconds'])
end

