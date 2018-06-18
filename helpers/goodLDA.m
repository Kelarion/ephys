function W = goodLDA(X,c)
% W = goodLDA(X,c)
%
% W is the set of eigenvectors of (Sw^-1)Sb associated with top C-1
% eigenvalues.
% These are the best linear discriminants of the clusters 'c' in X.
%
% X must nObs x nDim, and c is an nObs x 1 vector of class labels

c = c(:);
nGrp = length(unique(c));
nDim = size(X,2);
all_grp = unique(c);

bigmu = mean(X,1);

% get group means
tmpC = repmat(c,nDim,1) + repelem(nGrp*([1:nDim]' - 1),length(c));
mu_grp = clusterAverage(tmpC,X(:));
mu_grp = reshape(mu_grp,nGrp,nDim);

Sb = ((mu_grp - bigmu)'*(mu_grp - bigmu))/(nGrp - 1); % between-class covariance

Sw = zeros(nDim,nDim,nGrp); % within-class covariance
for iGrp = 1:nGrp
    Sw(:,:,iGrp) = cov(X(c == all_grp(iGrp),:));
end % should we average?
Sw = nansum(Sw,3);

[V,L] = eig(Sw\Sb);
[~,ord] = sort(diag(L),'descend');
notThese = ~real(diag(L)); useThese = ~notThese(ord); % only non-zero eigs

W = V(:,ord(useThese));
if size(W,2) >= nGrp
    W = W(:,1:(nGrp-1)); % trim off lesser discriminants
end
W = W';

end

%% to remove dependencies
function clusterQuantity = clusterAverage(clu, spikeQuantity)
% function clusterQuantity = clusterAverage(clu, spikeQuantity)
%
% get the average of some quantity across spikes in each cluster, given the
% quantity for each spike
%
% e.g.
% > clusterDepths = clusterAverage(clu, spikeDepths);
%
% clu and spikeQuantity must be vector, same size

% using a super-tricky algorithm for this - when you make a sparse
% array, the values of any duplicate indices are added. So this is the
% fastest way I know to make the sum of the entries of spikeQuantity for each of
% the unique entries of clu
[~, spikeCounts] = countUnique(clu);

% convert clu to indices, i.e. just values between 1 and nClusters.
[~,~,cluInds] = unique(clu);

% summation
q = full(sparse(cluInds, ones(size(clu)), double(spikeQuantity)));

% had sums, so dividing by spike counts gives the mean depth of each cluster
clusterQuantity = q./spikeCounts;

    function [values,instances] = countUnique(x)
        % function [values,instances] = countUnique(x)
        
        y = sort(x(:));
        p = find([true;diff(y)~=0;true]);
        values = y(p(1:end-1));
        instances = diff(p);
        
    end

end