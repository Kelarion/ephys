function shuffledArray = chopnotslop(X,P,dims)
% shuffledArray = chopnotslop(X,P,dims)
%
% Works on arbitrary (I think) N-dimensional arrays. 

% P is a matrix of size: size(X,dims(2))-by-size(X,dims(1)). 
% X will be shuffled along dims(2) independently for each of dims(1). 
% To do this, we rotate X around so that dims(1) is the columns of X, 
% and dims(2) is the rows; all others are trailing dimensions.

trailers = setdiff(find(size(X)>1),dims); 
X_twisted = permute(X,[fliplr(dims),trailers]);
nonsingles = find(size(X_twisted)>1);
dimLength = size(X_twisted);

assert(all(dimLength(dims) == size(P)),'P must be the same size as X in the plane being shuffled')

% The way we do it efficiently is by implicitly vectorizing X.
% Since each entry of P is the index of dims(2) along some other dimension, 
% we need to add to the corrected index the distance along vectorized X to 
% get to that other dimension. The following way is kind of werid, but it's
% so that this can generalize to any N-d array.
size_delta = dimLength;
size_delta([1 2]) = 1;
delta = repmat([1:dimLength(1)]',size_delta);

for iDim = trailers % now we need to shift each trailing dimension according
    inds = {};      % to how many elements would precede it in X(:)
    otherDims = setdiff(trailers,iDim);
    inds{1} = ':'; % weird matlab thing that works somehow
    inds{2} = 1; 
    inds(otherDims) = {':'};
    dd = prod(dimLength(nonsingles < iDim));  
    for ii = 2:dimLength(iDim)  % this will be added to all dimensions
        inds{iDim} = ii;  % and the next level of shifting is added later
        delta(inds{:}) = delta(inds{:}) + dd*[(ii-1)];
    end
end

shuffledArray = X_twisted(bsxfun(@plus,delta,(P-1)*dimLength(1)));
shuffledArray = permute(shuffledArray,[fliplr(dims),trailers]); % back to size(X)
