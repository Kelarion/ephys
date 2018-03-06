function [opto] = cleanAux(raux)
% Used to clean up contamination in the aux channel (due to some weird
% electrical issue we had once). Takes in just the auxillary channel, and
% will pull out the opto time stamps.



% first find out when the video is on, and restrict search to that time
[all_pks, all_locs] = findpeaks(double(raux)); 
first_peak = find(diff(all_pks) > 10^3, 1, 'first'); % heuristic
last_peak = find(diff(all_pks) > 10^3, 1, 'last');
vid_on = all_locs(first_peak):all_locs(last_peak);
pks = all_pks(diff(all_pks) > 10^3); 
locs = all_locs(diff(all_pks) > 10^3);

% then do the dumbest thing imagineable
nbin = median(diff(locs(1:100)));
M = movmean(raux,nbin);
opto = M > max(M)/3; % threshold the moving average


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