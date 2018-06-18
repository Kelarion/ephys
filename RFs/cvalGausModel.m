function [cvscore, modPars] = cvalGausModel(mdl,st,stimTimes,stimPos,params)
% [cvscore,modPars] = cvalGausModel(mdl,st,stimTimes,stimPos,params)
%
% k-fold cross-validation of RF model, with optional stochastic fitting.
%
% Inputs:
%   - mdl: struct with fields:
%       > String: one of {'Gaussian','DoG','DoG Unequal','Gabor'}
%           or {'flat','rf'} for null models
%       > func: function handle to the corresponding gaussian function
%       > fitIpsi: logical, if true then fits a second (gaussian) to the
%       ipsilateral hemifield, and adds this to the main fit (takes twice
%       as long, since it's fitting 2 models on each fold).
%       > ipsitype: 'same' if the ipsilateral model is the same as the
%       contralateral, or 'gaussian' if it should be a (unimodal) gaussian.
%   - st: spike times of neuron
%   - stimTimes: time of each dot/square appearing in the screen
%   - stimPos [ndot-by-2]: position of each dot/square given in (y,x)
%   - params: struct with fields:
%       > nFold: fold of CV (default 6)
%   	> nStoch: how many times to to run the fit, perturbing initial
%       conditions on each run (default 1, no perturbations)
%       > method: how to compare the model and spike train, can be:
%           - 'direct': (default) convert both to spike counts/spike
%           probability, rescale and take the R^2
%           - 'itskov': compute the Q metric (Itskov et al., 2008) on
%           predicted lambda and the spike train
%           - 'ZOH': replace fitted timecourse with zero-order hold
%           betweeen stimuli, and compute Q on this model (currently not
%           supported, will just use itskov).
%
% Output:
%   - cvscore: score for each fold, either Q or R^2 depending on method

snrfparams.useSVD = true;
snrfparams.makePlots = false;
snrfparams.fit2Dgauss = false;

if nargin < 5, params = struct; end
if ~isfield(params,'nFold'), params.nFold = 6; end
if ~isfield(params,'nStoch'), params.nStoch = 1; end
if ~isfield(params,'method'), params.method = 'direct'; end

bs_spk = 0.035;% time for spike binning, use median inter-stimulus time

if ~isfield(mdl,'fitIpsi'); mdl.fitIpsi = false; end
if ~isfield(mdl,'ipsitype'); mdl.ipsitype = 'gaussian'; end

yPos = unique(stimPos(:,1)); nY = length(yPos);
xPos = unique(stimPos(:,2)); nX = length(xPos);
% strange way to get all the possible (x,y) coords in the order I want
[X,Y] = meshgrid(xPos,yPos);
allPos = fliplr(reshape(permute(cat(3,X,Y),[3 1 2]),2,nX*nY)');

[t,ord] = sort(stimTimes);
p = stimPos(ord,:);
whenStim = unique(t);
nStim = length(whenStim);

% fit the main model first if what we want to cross-validate is the
% ipsilateral model
if mdl.fitIpsi
    roof = sparseNoiseRF_MA(st,t,p,snrfparams);
    bnds = setBounds(mdl,roof,xPos,yPos,1.5);
    [~, fit_mod] = stochfit(mdl,roof,xPos,yPos,bnds,20);
end

d = floor(nStim/params.nFold);
cvscore = zeros(params.nFold,1);
modPars = [];
for k = 1:params.nFold
    theseT = saturate([1:d] + d*(k-1),nStim);
    testStims = ismember(t,whenStim(theseT));
    fitStims = ~testStims;
    test_lims = [min(t(testStims)) max(t(testStims))];
    T = abs(diff(test_lims));
    
    % fit on part
    [RF,sts] = sparseNoiseRF_MA(st,t(fitStims),p(fitStims,:),snrfparams);
    
    if mdl.fitIpsi % support for an ipsilateral model
        ipsmdl = mdl;
        switch mdl.ipsitype
            case 'same'
                ipsmdl.String = mdl.String;
                ipsmdl.func = mdl.func;
            case 'gaussian'
                ipsmdl.String = 'gaussian';
                ipsmdl.func = @D2GaussFunctionRot;
        end
        ibnds = setBounds(ipsmdl,RF,-xPos,yPos,1.5); % xPos is flipped (b/c other hemifield)
        [x_fit, ipsi_mod] = stochfit(ipsmdl,RF-fit_mod,xPos,yPos,ibnds,params.nStoch);
        full_mod = fit_mod+ipsi_mod;
    else
        bnds = setBounds(mdl,RF,xPos,yPos,1.5);
        [x_fit, full_mod] = stochfit(mdl,RF,xPos,yPos,bnds,params.nStoch);
    end
    modPars = [modPars, x_fit(:)];
    
    % test on remaining
    [r_pred,t_pred] = linFiltStimulus(full_mod(:),sts,t(testStims),p(testStims,:),allPos);
    r_pred = rectify(r_pred + sts.bsl); % prediction is in units of deviation from baseline
    
    switch lower(params.method)
        case 'direct' % convert both to same
            [r_true,bins] = timestampsToBinned(st,0,bs_spk,test_lims);
            
            % interpolate prediction to same timeframe as actual rate
            pred = interp1(t_pred,r_pred,bins);
            
            actual = r_true;
            pred = pred*(mean(actual)/mean(pred)); % scaling fudge factor
            
            cvscore(k) = 1 - (std(actual-pred)^2)/(std(actual)^2); % R^2
        case 'itskov'
            % prediction quality 'Q', a la Itskov et al. (2008)
            actualSpikes = st(logical(WithinRanges(st,test_lims)));
            if ~isempty(actualSpikes)
%                 sp_inds = findNearestPoint(actualSpikes,t_pred);
                sp_inds = closest(t_pred(:),actualSpikes(:));
            else
                sp_inds = [];
            end
            bs = mean(diff(t_pred));
            cvscore(k) = (2*sum(r_pred(sp_inds)) - trapint(r_pred.^2,bs))/T;
        case 'zoh' % not supported yet
            warning('Sorry, ZOH this isn''t supported yet, just doing regular Itskov')
            actualSpikes = st(logical(WithinRanges(st,test_lims)));
            if ~isempty(actualSpikes)
%                 sp_inds = findNearestPoint(actualSpikes,t_pred);
                sp_inds = closest(t_pred(:),actualSpikes(:));
            else
                sp_inds = [];
            end
            bs = mean(diff(t_pred));
            cvscore(k) = (2*sum(r_pred(sp_inds)) - trapint(r_pred.^2,bs))/T;
            %             actualSpikes = st(logical(WithinRanges(st,test_lims)));
            %             sp_inds = discretize(actualSpikes,t_pred);
            %             % use Q again
            %             cvscore(k) = (2*sum(r_pred(sp_inds)) - norm(r_pred)^2);
        otherwise
            error('Invalid metric, I can''t compute that yet')
    end
    
end


end


%% numerical integrator
function X = trapint(x,bs)
% trapezoidal integration of x, with uniform spacing bs

X = (bs/2)*sum(x(1:end-1) + x(2:end));

end



