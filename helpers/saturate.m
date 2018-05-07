function p = saturate(r,z)
% Opposite of rectify function
%

p = r;
p(r>z) = z;
reshape(p,size(r));

end