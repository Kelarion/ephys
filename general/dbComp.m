function [inCommon, inds1, inds2] = dbComp(db1,db2)
% [inCommon, inds1, inds2] = dbComp(db1,db2)
% 
% Compare two db structs, return the common entries

foo1 = cat(1,{db1(:).mouse_name},{db1(:).date})';
foo2 = cat(1,{db2(:).mouse_name},{db2(:).date})';

inBoth = logical(prod(ismember(foo1,foo2),2));
inCommon = db1(inBoth);

inds1 = find(logical(prod(ismember(foo1,foo2),2)));
inds2 = find(logical(prod(ismember(foo2,foo1),2)));

end