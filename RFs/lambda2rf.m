function rf = lambda2rf(lamb,t,times,pos)


[~,~,whichPos] = unique(pos,'rows');
avg = eventLockedAvg(lamb,t,times,whichPos,[0 0.2]);
avg = permute(avg,[1 3 2]);

huh = unique(pos,'rows');
[~,ord] = sort(huh(:,2));

rf = reshape(max(avg(ord,:),[],2),10,36);