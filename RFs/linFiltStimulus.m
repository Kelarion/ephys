function [lambda,t_eval] = linFiltStimulus(RF,sts,stimTimes,stimPos,allPos)
% [lambda,t] = linFiltStimulus(RF,stimTimes,stimPos)
% 
% Model firing rate as a linear filtering of the stimulus. RF is the
% spatial component, and the temporal is in sts.timeCourse;

if ~exist('allPos','var')
    [X,Y] = meshgrid(unique(stimPos(:,2)),unique(stimPos(:,1)));
    nX = length(unique(stimPos(:,2))); nY = length(unique(stimPos(:,1)));
    allPos = fliplr(reshape(permute(cat(3,X,Y),[3 1 2]),2,nX*nY)');
end

[t,ord] = sort(stimTimes);
p = stimPos(ord,:);

[~,posID] = ismember(p,allPos,'rows');
[~,timID] = ismember(t,unique(t),'rows');
allStim = sparse(posID,timID,true, length(allPos),length(unique(t)));

lambda = (RF')*allStim;
% t_eval = unique(t);

% convolve the result with the timecourse
bs = mean(diff(sts.timeBins));
t_eval = (t(1)-bs*2):bs:(t(end)+max(sts.timeBins));
t_grps = findNearestPoint(unique(t),t_eval); % indices of evaluated points

lmb = sparse(ones(size(t_grps)),t_grps,lambda,1,length(t_eval));

ll_temp = sconv2(lmb,sts.timeCourse)*sts.scalar;
lambda = full(ll_temp(1:length(t_eval)));

end