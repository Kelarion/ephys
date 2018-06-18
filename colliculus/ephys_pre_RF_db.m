load('C:\Users\Matteo A\Documents\MATLAB\mine\sortRFfiles\bilats.mat')

clear db
shanks = fliplr(unique(bilats.whichCell(:,[3 2]),'rows'));
unique_k = unique(shanks(:,2));
nshank = length(unique_k);

for n = 1:nshank
    theseTags = shanks(shanks(:,2) == unique_k(n),1);
    db(n) = bilats.datasets(unique_k(n));
    db(n).tags = bilats.datasets(unique_k(n)).tags(theseTags);
    
    [expNums, ~, hasBlock] = ...
        dat.whichExpNums(db(n).mouse_name, db(n).date);
    
    db(n).noiseExp = expNums(find(hasBlock,1,'first'));
end

db = rmfield(db,{'tlExp','cwExp','passiveExp'});