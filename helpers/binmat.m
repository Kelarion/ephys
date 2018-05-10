function B = binmat(n,bs,noverlap)
% B = binmat(n,bs,noverlap)
%
% Make a matrix B such that B*x returns x summed in bins of size 'bs'. 
% B is size nbin x n, and x should be an n x 1 vector.

if ~exist('noverlap','var')
    noverlap = 0; 
elseif noverlap >= bs
    noverlap = 0; 
    warning('noverlap can''t be larger than bin size, making it 0');
end

nbin = ceil(n/(bs-noverlap));

cols = 1:n;
rws = repelem(1:nbin,bs);
rws = rws(1:n);

B = sparse(rws,cols,true,nbin,n); % need to add back support for overlap

% % bs = floor(n/nbin) + mod(n,nbin); % bin size
% B = zeros(nbin,n);
% for r = 1:nbin
%    i = (bs - noverlap)*(r-1)+ 1;
%    B(r,i:i+bs-1) = 1;
% end
% 
% B = B(:,1:n);
end