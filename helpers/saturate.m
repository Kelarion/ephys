function p = saturate(r,z)
% Opposite of rectify function
%
p = zeros(size(r));
for i = 1:length(r)
    if r(i) > z
        p(i) = z;
    else
        p(i) = r(i);
    end
end

end