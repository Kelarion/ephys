function [whichRows, stepSize] = findDiagBlocks(BD,k,algo)
% [rows] = findDiagBlocks(BD[,k,algo])
% 
% Matrix must be symmetric, gives rubbish if it isn't actually
% block-diagonal in structure

if nargin <3, algo = 'svd'; end

switch algo
    case 'eigs' 
        [V, D] = eigs((BD > mean(BD(:)))+0); 
    case 'svd' % works better, presumably because squaring the
        [V, D, ~] = svd(BD); % matrix amplifies block structure?
end

if nargin<2 || isempty(k)
    lambda = diag(D)/sum(diag(D));
    k = sum(lambda >= max(lambda)/2);
end

whichRows = zeros(k,1);
stepSize = zeros(k,1);
for iRow = 1:k
    smoothV = smooth(zscore(V(:,iRow)),0.1,'loess'); % cook heavily
    [~, rStep] = max(abs(diff(smoothV)));
    whichRows(iRow) = rStep + 1; % because we used 'diff'
    conf = (mean(V(rStep:end,iRow))-mean(V(1:rStep,iRow)))/range(V(:));
    stepSize(iRow) = conf;
end