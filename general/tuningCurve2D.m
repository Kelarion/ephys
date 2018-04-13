function [resp,errs,numEvs,whichEvs,psth,bins,rast] = tuningCurve2D(spikeTimes,eventTimes,var1,var2,win,binsize)
% [resp,errs,numEvs,whichEvs,psth,bins,rast] = tuningCurve2D(spikeTimes,eventTimes,var1,var2,win[,binsize])
%
% default bin size is 1 ms

if nargin<6;binsize = 0.001;end

[~,bins,~,~,~,ba] = psthAndBA(spikeTimes,eventTimes,win,binsize);

v1 = unique(var1);
v2 = unique(var2);
resp = zeros(length(v2),length(v1));
errs = zeros(length(v2),length(v1));
psth = zeros(length(v2)*length(v1),round(diff(win)/binsize));
numEvs = zeros(length(v2),length(v1));
for ii = 1:length(v2)
    for jj = 1:length(v1)
        theseEvents = var1 == v2(ii) & var2 == v1(jj);
        p = mean(ba(theseEvents,:)./binsize,1);
        
        numEvs(ii,jj) = sum(theseEvents);
        resp(ii,jj) = mean(sum(ba(theseEvents,:),2)/diff(win));
        errs(ii,jj) = std(sum(ba(theseEvents,:),2)/diff(win))./sqrt(sum(theseEvents));
        psth(ii + 4*(jj-1),:) = p;
    end
end

whichEvs = zeros(length(v2),length(v1),2);
[x,y] = meshgrid(v2,v1);
whichEvs(:,:,1) = x;
whichEvs(:,:,2) = y;

[tr,b] = find(ba);

[rasterX,yy] = rasterize(bins(b));
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3);
rasterScale = 6;
rasterY(2:3:end) = rasterY(2:3:end)+rasterScale;

rast = [rasterX(:) rasterY(:)];