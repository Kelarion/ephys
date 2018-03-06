function [eta, peri] = etrigav(signal, event, window)
% [eta peri] = etrigav(signal, event, window)
%

nSigs = min(size(signal));

evOn = find(diff(event) > 0); % get the onset and offset times of the event

if nargin < 3
    evOff = find(diff(event) < 0);% 1 x nEvs in this event
    longest_event = max(evOff - evOn);
    winSize = longest_event + 50;
    delay = 0;
else
    winSize = window(2)-window(1)+1;
    delay = abs(window(1));
end

eta = zeros(nSigs, winSize); 
peri = zeros(length(evOn), winSize, nSigs);
for iSig = 1:nSigs
    for jEv = 1:length(evOn)
        ev_start = rectify(evOn(jEv) - delay,1);
        ev_fin = ev_start + winSize - 1;
        if ev_fin > length(event)
            peri(end,:) = [];
            continue
        end
        peri(jEv, 1:winSize, iSig) = signal(iSig, ev_start:ev_fin);
    end
    av = nanmean(peri(:, :, iSig), 1);
    eta(iSig, 1:winSize) = av;
end

if nargout < 1
   dims = factor(nSigs); % dimensions of the subplot
   if length(dims) > 2
       h = ceil(length(dims)/2);
       nrow = prod(dims(1:h));
   else 
       nrow = dims(1);
   end
   for iSig = 1:nSigs
       subplot(nrow, nSigs/nrow, iSig)
       title(sprintf('Neuron %d', iSig))
       p1 = plot(peri(:,:,iSig)','Color',[0.6 0.6 0.6]);
       hold on
       p2 = plot(eta(iSig,:),'Color','k','LineWidth',4);
   end
   legend([p1(1) p2],'Individual trials', 'Average')
end

end

function p = rectify(r,z)

if ~exist('z','var')
   z = 0; 
end
if r < z 
    p = z;
else
    p = r;
end

end