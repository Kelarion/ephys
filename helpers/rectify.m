function p = rectify(r,z)

if ~exist('z','var')
   z = 0; 
end
p = zeros(size(r));
for i = 1:length(r)
    if r(i) < z
        p(i) = z;
    else
        p(i) = r(i);
    end
end

end