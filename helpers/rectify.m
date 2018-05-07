function p = rectify(r,z)
% p = rectify(r,z)
%
% Returns p = {r, if r >= z
%             {0, is r <  z

if ~exist('z','var')
   z = 0; 
end

p = r;
p(r<z) = z;
reshape(p,size(r));

end