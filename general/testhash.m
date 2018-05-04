%% params
n= 120;
nclu_max = 18:1100;
db_max = [4 6 8 10 12];
bases = 3:12;
wc2id = @(WC,b) WC*[b^2 b^1 b^0]'; % basis conversion
id2wc = @(id,b) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)];
%% generate fake serials 
fake = zeros(n,3,length(nclu_max),length(db_max));
for d = 1:length(db_max)
    for k = 1:length(nclu_max)
        ftags = [];
        dd = db_max(d);
        for ii = 1:dd
            ftags = [ftags; sort(randi(4,n/dd,1))];
        end
        fcids = zeros(size(ftags));
        for ii = 1:4
            nn = diff(find(abs(diff([0; ismember(ftags,ii); 0]))));
            nn = nn(1:2:end);
            fakecids = [];
            for jj = 1:length(nn)
                fakecids = [fakecids sort(randperm(nclu_max(k),nn(jj)))];
            end
            fcids(ismember(ftags,ii))  = fakecids;
        end
        fake(:,:,k,d) = [fcids ftags repelem([1:dd]',n/dd)];
    end
end
%% test
dmntly = zeros(length(bases),length(nclu_max),length(db_max));
invtblty = zeros(length(bases),length(nclu_max),length(db_max));
maxval = zeros(length(bases),length(nclu_max),length(db_max));
whenUnique = zeros(length(nclu_max),length(db_max));
whenInvertible = zeros(length(nclu_max),length(db_max));
for d = 1:length(db_max)
    for k = 1:length(nclu_max)
        for ii = 1:length(bases)
            foo = wc2id(fake(:,:,k),bases(ii));
            dmntly(ii,k,d) = length(unique(foo));
            ba = id2wc(foo,bases(ii));
            invtblty(ii,k,d) = sum(sum(ba == fake(:,:,k)));
            maxval(ii,k,d) = max(foo);
        end
        whenUnique(k,d) = bases(find(dmntly(:,k) == n,1,'first'));
        whenInvertible(k,d) = bases(find(invtblty(:,k) == n*3,1,'first'));
    end
end

%%
figure
subplot(1,2,1)
imagesc(nclu_max,db_max,whenUnique)
xlabel('size of first ''dimension''')
ylabel('size of last ''dimension''')
title('minimum for unique map')

subplot(1,2,2)
imagesc(nclu_max,db_max,whenInvertible)
xlabel('size of first ''dimension''')
title('minimum for invertible map')

colorbar


