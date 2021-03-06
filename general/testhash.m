%% params
nclu = 1000;
CID_max = [600:100:1200];
dset_max = [4:4:55];
tag_max = 4;
bases = 3:56;

wc2id = @(WC,b) WC*[b^2 b^1 b^0]'; % basis conversion
id2wc = @(id,b) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)];
%% generate fake serial IDs
% 'fake' is an nClu-by-3 matrix, generated for each combination of maximum
% CID value and maximum dset value (max tag value is assumed to be <= dset)
fake = zeros(nclu,3,length(CID_max),length(dset_max));
for d = 1:length(dset_max)
    for k = 1:length(CID_max)
        % simulate heirarchical structure of the ID matrices
        fake_dsets = sort(discretize(randperm(nclu),dset_max(d)));
        fake(:,3,k,d) = fake_dsets;
        for iset = 1:dset_max(d)
            thisSet = (fake_dsets == iset);
            fake(thisSet,2,k,d) = sort(discretize(randperm(sum(thisSet)),tag_max));
            for itag = 1:tag_max
                thisTag = (fake(:,2,k,d) == itag) & thisSet(:);
                fake(thisTag,1,k,d) = sort(randperm(nclu,sum(thisTag)))';
            end
        end
        
        %         % use random pairings
        %         fake(:,:,k,d) = [randi(CID_max(k),nclu,1) ...
        %             randi(dset_max(d),nclu,1) ...
        %             randi(tag_max,nclu,1)];
    end
end
%% test
% here we pass 'fake' through the mapping for each combination of maximum
% CID and dset value, for different choice of base, and find when that
% mapping is unique; then we do the inverse mapping and find when that
% gives us back the original 'fake' IDs.
dmntly = zeros(length(bases),length(CID_max),length(dset_max));
invtblty = zeros(length(bases),length(CID_max),length(dset_max));
maxval = zeros(length(bases),length(CID_max),length(dset_max));
whenUnique = zeros(length(CID_max),length(dset_max));
whenInvertible = zeros(length(CID_max),length(dset_max));
for d = 1:length(dset_max)
    for k = 1:length(CID_max)
        for ii = 1:length(bases)
            foo = wc2id(fake(:,:,k,d),bases(ii));
            dmntly(ii,k,d) = length(unique(foo));
            
            ba = id2wc(foo,bases(ii));
            invtblty(ii,k,d) = sum(any(ba ~= fake(:,:,k,d)));
            maxval(ii,k,d) = max(foo);
        end
        original = length(unique(fake(:,:,k,d),'rows'));
        whenUnique(k,d) = bases(find(dmntly(:,k,d) == original,1,'first'));
        whenInvertible(k,d) = bases(find(invtblty(:,k,d) == 0,1,'first'));
    end
end

%%
figure
subplot(1,2,1)
imagesc(dset_max,CID_max,whenUnique)
xlabel('max value of dset')
ylabel('max value of CID')
title('minimum base for unique map')

subplot(1,2,2)
imagesc(dset_max,CID_max,whenInvertible)
xlabel('max value of dset')
title('minimum base for an invertible map')

colorbar


