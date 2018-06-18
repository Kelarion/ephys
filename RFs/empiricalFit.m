function [score, pred, actual] = empiricalFit(RF,stats,snrf,st,metric)
% score = empiricalFit(RF,stats,snrf,spikeTimes,metric)
%
% run the linear prediction of firing rate given RF and timecourse, compute
% the performance relative to specified metric (default R^2).
%
% stats is the same as returned by sparseNoiseRF_MA, must have baseline,
% scalar and timecourse fields (bsl, scalar, and timeCourse respectively).
%
% snrf 

if nargin < 5, metric = 'direct'; end

bs_spk = 0.05;% median inter-stimulus time
test_lims = [min(st) max(st)];

[r_pred, t_pred] = linFiltStimulus(RF(:),stats,snrf.stimTimes_local,snrf.stimPosition);
r_pred = rectify(r_pred + stats.bsl);

switch metric
    case 'direct'
        [r_true,bins] = timestampsToBinned(st,0,bs_spk,test_lims);
%         r_pred = rectify(r_pred + mean(r_true));
        
        % interpolate prediction to same timeframe as actual rate
        pred = interp1(t_pred,r_pred,bins);
        
        actual = r_true;
        pred = pred*(mean(actual)/mean(pred)); % rescale to have same means
        score = 1 - (std(actual-pred)^2)/(std(actual)^2);
    case 'itskov'
        r_pred = rectify(r_pred + stats.bsl);
        if ~isempty(st)
            %                 sp_inds = findNearestPoint(actualSpikes,t_pred);
            sp_inds = closest(t_pred(:),st(:));
        else
            sp_inds = [];
        end
        bs = mean(diff(t_pred));
        score = (2*sum(r_pred(sp_inds)) - trapint(r_pred.^2,bs))/range(st);
        pred = r_pred;
        actual = [];
    otherwise
        error('Invalid metric, I can''t compute that yet')
end


end

%% numerical integrator
function X = trapint(x,bs)
% trapezoidal integration of x, with uniform spacing bs

X = (bs/2)*sum(x(1:end-1) + x(2:end));

end
