function [opto] = cleanAux(raux)
% Used to clean up contamination in the aux channel (due to some weird
% electrical issue we had once). Takes in just the auxillary channel, and
% will pull out the opto time stamps (in terms of indices of input).

t = 1:length(raux);
% first find out when the video is on, and restrict search to that time
[all_pks, all_locs] = findpeaks(double(raux)); 
locs = all_locs(diff(all_pks) > 10^3); % look only when video is on

% then do the dumbest thing imagineable
nbin = median(diff(locs(1:100)));
M = movmean(raux,nbin);
isopto = M > max(M)/3; % threshold the moving average
optoon = find(diff([0 isopto])>0);
optooff = find(diff([0 isopto])<0);

opto = [optoon(:) optooff(:)];

tooSoon = find(diff(opto,[],2) < 50); % events occuring too soon after the last
opto(tooSoon-1,2) = opto(tooSoon,2); % merge with adjacent events
opto(tooSoon,:) = [];
opto = t(opto);

% % experimental
% if nargout > 2 
%     dopto = diff(opto); % still 1 x nPks
%     opton = find(dopto > 0);
%     optoff = find(dopto < 0);
%     for iOpt = 1:sum(opton)
%         before = locs(opton(iOpt)) - 3;
%         after = locs(optoff(iOpt)) + 1;
%         predicted_vid = repmat();
%         newraux(before:after)
%     end
% end

end