% load('C:\Users\Matteo A\Documents\MATLAB\mine\ephys\cerebellum\allrips.mat')
load('C:\Users\Matteo A\Documents\MATLAB\mine\ephys\cerebellum\gt_easy.mat')
load('C:\Users\Matteo A\Documents\MATLAB\mine\ephys\cerebellum\gt_hard.mat')

wc2id = @(WC) WC*[5^2 5^1 5^0]';
id2wc = @(id) [floor(id/(5^2)) mod(floor(id/5),5) mod(id,5)];

dset1 = gt_hard.dset; % hard dataset
t_max1 = gt_hard.t_max;

dset2 = gt_easy.dset; % easy dataset
t_max2 = gt_easy.t_max;

db = R.datasets;
wc = [dset1;dset2];

oscID = wc2id(R.osc.whichCell);
cgs = unique(oscID);
nClu = length(cgs);
%% get datasets
for iSet = 1:2
    k = wc(iSet,3);
    d = wc(iSet,2);
    parentDir = [ksRoot '\' db(k).name '\'];
    subDir = dir([parentDir '*' db(k).date]);
    
    if length(subDir) > 1
        subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
        if isempty(subDir)
            disp('can''t find it ¯\_(?)_/¯')
        end
    end
    dataDir = [parentDir subDir.name '\'];
    
    sp(iSet) = loadJRCdir(dataDir);
end

algo1 = R.osc.Times(oscID == wc2id(dset1),:);
algo1 = algo1(algo1(:,2) < t_max1,:);
algo2 = R.osc.Times(oscID == wc2id(dset2),:);
algo2 = algo2(algo2(:,2) < t_max2,:); 

%% establish ground truth
gt1 = gt_hard.Times;
gt2 = gt_easy.Times;

%% hard dataset

% first, which ground truth events did it correctly find
sameEv = [];
for iEv = 1:size(algo1,1)
    overlaps = (algo1(iEv,1)<=gt1(:,2)) & (algo1(iEv,2)>=gt1(:,1));
    if any(overlaps)
        sameEv = [sameEv; find(overlaps)]; % indices of GT events that were found
    end
end

gt_found = ismember(1:size(gt1,1),unique(sameEv)); % hits
misses = ~gt_found; % detected events that were not real

% next, which detected events did not correspond to GT event
sameEv = [];
for iEv = 1:size(gt1,1)
    overlaps = (gt1(iEv,1)<=algo1(:,2)) & (gt1(iEv,2)>=algo1(:,1));
    if any(overlaps)
        sameEv = [sameEv; find(overlaps)]; % indices of detected events that were real
    end
end

algo_right = ismember(1:size(algo1,1),unique(sameEv));
FA = ~algo_right;
hits = min([sum(algo_right) sum(gt_found)]); % in case of double-detection

if (sum(algo_right) > sum(gt_found))
    inds = findNearestPoint(gt1(gt_found,1),algo1(algo_right,1));
    foo = algo1(algo_right,1);
    hard_onsets = [foo(inds) gt1(gt_found,1)];
elseif (sum(algo_right) < sum(gt_found))
    inds = findNearestPoint(algo1(algo_right,1),gt1(gt_found,1));
    foo = gt1(gt_found,1);
    hard_onsets = [algo1(algo_right,1) foo(inds)];
else
    hard_onsets = [algo1(algo_right,1) gt1(gt_found,1)];
end
    
hard_perf = [hits sum(FA);sum(misses) NaN];
hard_missed = gt1(misses,:);
hard_FA = algo1(FA,:);

%% easy dataset

% first, which ground truth events did it correctly find
sameEv = [];
for iEv = 1:size(algo2,1)
    overlaps = (algo2(iEv,1)<=gt2(:,2)) & (algo2(iEv,2)>=gt2(:,1));
    if any(overlaps)
        sameEv = [sameEv; find(overlaps)]; % indices of GT events that were found
    end
end

gt_found = ismember(1:size(gt2,1),unique(sameEv)); % hits
misses = ~gt_found; % detected events that were not real

% next, which detected events did not correspond to GT event
sameEv = [];
for iEv = 1:size(gt2,1)
    overlaps = (gt2(iEv,1)<=algo2(:,2)) & (gt2(iEv,2)>=algo2(:,1));
    if any(overlaps)
        sameEv = [sameEv; find(overlaps)]; % indices of detected events that were real
    end
end

algo_right = ismember(1:size(algo2,1),unique(sameEv));
FA = ~algo_right;
hits = min([sum(algo_right) sum(gt_found)]); % in case of double-detection

if (sum(algo_right) > sum(gt_found)) % onset times of true-positive events
    inds = findNearestPoint(gt2(gt_found,1),algo2(algo_right,1));
    foo = algo2(algo_right,1);
    ez_onsets = [foo(inds) gt2(gt_found,1)];
elseif (sum(algo_right) < sum(gt_found))
    inds = findNearestPoint(algo2(algo_right,1),gt2(gt_found,1));
    foo = gt2(gt_found,1);
    ez_onsets = [algo2(algo_right,1) foo(inds)];
else
    inds = findNearestPoint(algo2(algo_right,1),gt2(gt_found,1));
    ez_onsets = [algo2(algo_right,1) gt2(gt_found,1)];
end

ez_perf = [hits sum(FA);sum(misses) NaN];
ez_missed = gt2(misses,:);
ez_FA = algo2(FA,:);

%% summary plots
% total performance
figure('units','normalized','position',[0.0297 0.4741 0.1719 0.2463])
cumperf = ez_perf + hard_perf;
imagesc(cumperf)
yticks([1 2]);xticks([1 2])
yticklabels({'Y','N'})
xticklabels({'Y','N'})
hold on
text([1 1 2 2],[1 2 1 2],split(num2str(cumperf(:)')),'fontweight','bold','fontsize',14,...
    'horizontalalignment','center')
xlabel('Event exists')
ylabel('Event detected')
title('Overall performance')

%% which sorts of events were detected/missed/false-positive
allMissed = {hard_missed,ez_missed};
allFound = {hard_onsets(:,2),ez_onsets(:,2)};
AVG_found = {};
peri_found = {};
AVG_missed = {};
peri_missed = {};
for iSet = 1:2
    k = wc(iSet,3);
    d = wc(iSet,2);
    parentDir = [ksRoot '\' db(k).name '\'];
    subDir = dir([parentDir '*' db(k).date]);
    
    if length(subDir) > 1
        subDir = subDir(contains({subDir(:).name},db(k).depth{d}));
    end
    dataDir = [parentDir subDir.name '\'];
    
    cid = wc(iSet,1);
    rawDat = memmapfile([dataDir sp(iSet).dat_path], 'Format',  ...
        {sp(iSet).dtype, [sp(iSet).n_channels_dat sp(iSet).nSampDat], 'x'});
    
    chan = sp(iSet).chanMap(sp(iSet).mainChannel(sp(iSet).cids == cid));
    lief = getLFP(rawDat.Data.x(chan,:),sp(iSet).sample_rate,filt_low);
    tdat = [1:size(rawDat.Data.x,2)]/sp(iSet).sample_rate;
    
    lckdTimes = phaseLock(allFound{iSet},lief,-pi,sp(iSet).sample_rate);
    [ota,~,peri] = eventLockedAvg(double(rawDat.Data.x(chan,:)), ...
        tdat,lckdTimes,ones(size(lckdTimes)),[0 0.5]);
    AVG_found{iSet} = permute(ota,[3 2 1]) - median(ota,3); % because the output is sill
    peri_found{iSet} = permute(peri,[3 1 2]) - median(peri,3)';
    
    lckdTimes = phaseLock(allMissed{iSet},lief,-pi,sp(iSet).sample_rate);
    [ota,t,peri] = eventLockedAvg(double(rawDat.Data.x(chan,:)), ...
        tdat,lckdTimes,ones(size(lckdTimes)),[0 0.5]);
    AVG_missed{iSet} = permute(ota,[3 2 1]) - median(ota,3); % because the output is sill
    peri_missed{iSet} = permute(peri,[3 1 2]) - median(peri,3)';
end



