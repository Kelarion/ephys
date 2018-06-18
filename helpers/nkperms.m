function [P,pie,antipie]= nkperms(V,n,k,isUnique)
% [P,pi,pi_inv] = nkperms(V,n,k,isUnique)
%
% It's like 'perms' but you can specify how many permutations (n) and the
% degree of those permutations (k: k <= length(V)). So, the result will be
% be an n-by-k matrix (as opposed to length(V)!-by-k, for all the
% possible permutations). 
% 
% Inputs:
%   V:      vector
%   n:      integer
%   k:      integer <= length(V).
%
% Outputs:
%   P:      n-by-k, rows are a random pi(V)
%   pi:     the permutations that take V -> P, i.e. V(pi) -> P
%   pi_inv: the inverse mapping, i.e. P(i,pi_inv(i,:)) -> V for all i. To 
%           efficiently implement this, try chopnotslop
%               > NB, this is only defined when k == length(V), otherwise 
%               it spits out a NaN array. 
%
% To do:
%   add input isUnique: if true, then each permutation is unique and n must 
%           be <= length(V)!.(default is false)

if nargin<4, isUnique = false; end
if nargin<3, k = length(V); end

pie = zeros(n,k);
antipie = nan(n,k);
for nn = 1:n
    if ~isUnique
        pie(nn,:) = randperm(length(V),k);
        if k == length(V)
            % 'a' is the order of the elements of V, 'b' is how those
            % elements are ordered in P. 
            [~,a] = sort(V);
            [~,b] = sort(V(pie(nn,:)));
            antipie(nn,:) = b(a);
        end
    else 
        %add this later
    end
end

P = V(pie);