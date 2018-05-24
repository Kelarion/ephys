%  C = closest(A,B)
%  Inputs:  A = full real double vector (reference), sorted ascending
%           B = full real double vector (test),      sorted ascending
%  Output:  C = index of closest value in A to value in B
%           i.e., A(C(i)) is the closest value in A to B(i)
%           Ties are resolved in favor of higher index
%           Is not currently coded to handle Inf's and NaN's consistently
%           C is the same size as B
%
% Programmer:  James Tursa