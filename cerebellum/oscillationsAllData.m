%% load 
load('C:\Users\Matteo\Documents\MATLAB\ephys\cerebellum\allrips.mat')

% To convert from the 3-number identification (which correspond to cell ID,
% tag, and dataset) to a single unique cell ID, I use the below mappings. 
% They're equivalent to a conversion between base 5 and 10. For some 
% reason, using base-5 makes this mapping always invertible. (I checked 
% this with randomly-generated ID matrices where the 'digits' could vary 
% between 1 and 2000.)
wc2id = @(WC,b) WC*[b^2 b^1 b^0]';
id2wc = @(id,b) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)];

CID = wc2id(R.whichCell,5);
cgs = unique(CID);
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

mags = sqrt(optoMag.*spontMag); % determine the magnitude of each point

dtheta = spontPhs - optoPhs;
% dtheta = min(abs([dtheta dtheta+2*pi dtheta-2*pi]),[],2);
phsDiff = mags.*exp(1i*dtheta);
phsors = reshape([zeros(size(phsDiff)) phsDiff nan(size(phsDiff))].',length(inboth)*3,1);

% plot
figure
subplot(2,1,1)
scatlab(optoPhs,spontPhs,inboth,48*mags,spontFreq,'filled')
samelims
xlabel('phase during stimulation')
ylabel('spontaneous phase')
grid on
hold on;plot(xlim,xlim)
plot(xlim,xlim - 2*pi)
plot(xlim,xlim + 2*pi)

subplot(2,1,2)
scatlab(real(phsDiff),imag(phsDiff),inboth,[],spontFreq,'filled')
samelims
set(gca,'XLim',[-1, 1])
set(gca,'YLim',[-1, 1])
hold on
plot(real(phsors),imag(phsors))
grid on

c = colorbar;
c.Label.String = 'mean spontaneous frequency (Hz)';

% %% a plot
% figure; 
% cm = jet(length(unique(spontFreq)));
% for ii = 1:length(inboth)
%     arc = spontMag(ii).*exp(1i*linspace(spontPhs(ii),spontPhs(ii)+dtheta(ii),100));
%     plot(real(arc),imag(arc),'linewidth',2*spontMag(ii))
%     hold on
%     scatter(real(exp(1i*spontPhs(ii))),imag(exp(1i*spontPhs(ii))),'filled')
%     scatter(real(exp(1i*optoPhs(ii))),imag(exp(1i*optoPhs(ii))))
% end

%% per-event analysis 



%% videos
d = 2;
k = 1;
deez = ~R.isOpto & R.whichCell(:,2) == d & R.whichCell(:,3) == k;
cams = 0;

parentDir = [ksRoot '\' db(k).name '\'];
subDir = dir([parentDir '*' db(k).date]);

if length(subDir) > 1
    subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
    if isempty(subDir)
        disp('can''t find it ¯\_(?)_/¯')
    end
end
dataDir = [parentDir subDir.name '\'];

sp = loadJRCdir(dataDir,true);
ft = sp.frameTimes/1000;
offset = ft(1);

et = unique(R.interval(deez,:),'rows');
tooSoon = find(diff([0;et(:,1)]) < 0.1);
while ~isempty(tooSoon) % events occuring too soon after the last
    et(tooSoon-1,:) = (et(tooSoon-1,:) + et(tooSoon,:))/2;
    et(tooSoon,:) = [];  % merge with adjacent events
    tooSoon = find(diff([0;et(:,1)]) < 0.01);
end

iCam = cams(1);
whichCam = sprintf('bias_video_cam_%d_',iCam);
vidDir = dir([dataDir whichCam '*.avi']); % get the video file
vidPlace = [vidDir.folder '\' vidDir.name];

vidObj = VideoReader(vidPlace);

%% play
win = [-0.5 0.5];

v = VideoWriter('test.avi');
v.FrameRate = vidObj.FrameRate/5;
open(v)
tic
for iEv = 1:size(et,1)
    t_start = et(iEv,1)-offset;
    t_fin = et(iEv,2)-offset;
    if t_start<0, continue;end
    vidObj.CurrentTime = rectify(t_start+win(1));

    figure('visible','off'); % 'units','normalized','position',[0.1240 0.2889 0.6995 0.5250],
    while vidObj.CurrentTime < t_fin+win(2)
        vidFrame = readFrame(vidObj);
        image(vidFrame);
        hold on
        text(vidObj.width,vidObj.height,num2str(vidObj.CurrentTime),...
            'color','yellow','horizontalalignment','right', ...
            'verticalalignment','bottom','fontsize',15)
        text(0,0,['Event #: ' num2str(iEv)],'color','yellow',...
            'horizontalalignment','left','verticalalignment','top',...
            'fontsize',13)
        if (vidObj.CurrentTime > t_start) && (vidObj.CurrentTime < t_fin)
            scatter(400,50,800,'y','filled')
        end
        hold off
        set(gca,'Visible','off');
        frem = getframe(gca);
        writeVideo(v,frem);
        %     pause(1/vidObj.FrameRate);
    end
end
close(v)
toc
