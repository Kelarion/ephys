function [x,fitRF] = stochfit(mdl,rf,xPos,yPos,bounds,nStoch)
% bounds must be a structure with fields x0,lb,ub,A_ineq,b_ineq and nlin
% (this is returned from setBounds) to be fed as the constraints of fmincon
% will return NaN for x if model is 'flat' or 'rf',

switch lower(mdl.String)
    case 'flat'
        x = NaN;
        fitRF = zeros(size(rf));
        return
    case 'rf'
        x = NaN;
        fitRF = rf;
        return
end

opts = optimset('Display','off','Algorithm','sqp');% is this the best one to use?

nX = length(xPos); nY = length(yPos);
[X,Y] = meshgrid(xPos,yPos);
xdat = cat(3,X,Y);

% errFun = @(x)norm(reshape(mdl.func(x,xdat),nX*nY,1) - rf(:),1); % try L1
errFun = @(x)norm(reshape(mdl.func(x,xdat),nX*nY,1) - rf(:))^2;
[x0,lb,ub,A_ineq,b_ineq,nlincoln] = readBounds(bounds);

err = zeros(1,nStoch);
mods = zeros(nY,nX,nStoch);
prms = zeros(length(x0),nStoch);
inits = zeros(length(x0),nStoch);
for iInit = 1:nStoch
    if iInit == 1
        eps = 0;
    else 
        eps = perturb(mdl,bounds);
    end
    
    new_inits = x0 + eps;
    tmpx = fmincon(errFun,new_inits,A_ineq,b_ineq,[],[],lb,ub,nlincoln,opts);
    tmpfit = mdl.func(tmpx,xdat);
    
    err(iInit) =  errFun(tmpx);
    mods(:,:,iInit) = tmpfit;
    prms(:,iInit) = tmpx;
    inits(:,iInit) = new_inits;
end
[~,best_model_for_now] = min(err);
fitRF = mods(:,:,best_model_for_now);
x = prms(:,best_model_for_now);

end


%% perturbations for stochastic fitting
function epsilon = perturb(mdl,bounds)
% because we primarily want ot perturb the means

opts = optimset('Display','off');
[x0,lb,ub,A_ineq,b_ineq] = readBounds(bounds);

f = -ones(size(x0)); % find the maximum possible perturbation which satisfies constraints
[eps_max,~,isFeas] = linprog(f,A_ineq,(b_ineq-A_ineq*x0(:)),[],[],lb(:) - x0(:),ub(:) - x0(:),opts);
if isFeas>0          % i.e. find argmax_sum(eps): A_ineq*(x0+eps) <= b_ineq
    eps_max = eps_max(:)';
else
    eps_max = zeros(size(x0));
end

switch lower(mdl.String) % model-specific discounting of perturbations
    case 'gaussian'      % larger -> less variation
        dis = [100, 4,  4,  6, 10, 5];
        w   = [1    7   5   7   5  pi/12]; % units of perturbation 
    case 'dog'
        dis = [100, 17, 10, 6, 10, 100, 100, 100, 5];
        w   = [1    7   5   7   5  1    5    5   pi/12];
    case 'doug'
        dis = [100, 17, 10, 6, 10, 100, 30, 30,7,7, 5];
        w   = [1    7   5   7   5  1    5    5 7 7  pi/12];
    case 'gabor'
        dis = [100, 17, 10, 6, 10, 200, 30, 5];
        w   = [1    7   5   7   5  pi/36  pi/36 pi/12];
end

epsilon = (w./dis).*randn(size(w)); % perturb by some proportion of w 
violators = A_ineq*(x0(:)+epsilon(:)) > b_ineq; % see which ones went too far
tooHigh = ((x0+epsilon) > ub);
tooLow = (x0+epsilon < lb);

epsilon(violators) = eps_max(violators);
epsilon(tooHigh) = ub(tooHigh) - x0(tooHigh);
epsilon(tooLow) = lb(tooLow) - x0(tooLow);

end 

%% 
function [x0,lb,ub,A_ineq,b_ineq,nlincon] = readBounds(bounds)

x0 = bounds.x0;
lb = bounds.lb;
ub = bounds.ub;
A_ineq = [];
b_ineq = [];
nlincon = [];
if isfield(bounds,'A_ineq'), A_ineq = bounds.A_ineq; end
if isfield(bounds,'b_ineq'), b_ineq = bounds.b_ineq; end
if isfield(bounds,'nlin'), nlincon = bounds.nlin; end

end


