%% load 
load('C:\Users\Matteo\Documents\UNI\University of College\Ken class\allrips.mat')
R.isOpto = logical(R.isOpto);

% To convert from the 3-number identification (which correspond to cell ID,
% tag, and dataset) to a single unique cell ID, I use the below mappings. 
% They're equivalent to a conversion between base 5 and 10. For some 
% reason, using base-5 makes this mapping always invertible. (I checked 
% this with randomly-generated ID matrices where the 'digits' could vary 
% between 1 and 2000.)
wc2id = @(WC,b) WC*[b^2 b^1 b^0]';
id2wc = @(id,b) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)];

CID = wc2id(R.whichCell,5);

%% mean phase comparisons
[~,~,cidOp] = unique(CID(R.isOpto));
[optoPhs,optoMag] = splitapply(@circmean,R.spkPhase(R.isOpto),cidOp);

[~,~,cidSp] = unique(CID(~R.isOpto));
[spontPhs,spontMag] = splitapply(@circmean,R.spkPhase(~R.isOpto),cidSp);
spontFreq = splitapply(@mean,R.freq(~R.isOpto),cidSp);

cidinds = 1:length(unique(CID));
inboth = cidinds(ismember(cidinds,intersect(cidSp,cidOp)));
spontPhs = spontPhs(inboth);
spontMag = spontMag(inboth);
spontFreq = spontFreq(inboth);
optoPhs = optoPhs(inboth);
optoMag = optoMag(inboth);

mags = sqrt(optoMag.*spontMag);

dtheta = spontPhs(inboth) - optoPhs(inboth);
% dtheta = min(abs([dtheta dtheta+2*pi dtheta-2*pi]),[],2);
phsDiff = mags.*exp(1i*dtheta);
phsors = reshape([zeros(size(phsDiff)) phsDiff nan(size(phsDiff))].',length(inboth)*3,1);

% plot
figure
subplot(2,1,1)
scatter(optoPhs,spontPhs,48*mags,spontFreq,'filled')
samelims
xlabel('phase during stimulation')
ylabel('spontaneous phase')
grid on
hold on;plot(xlim,xlim)
plot(xlim,xlim - 2*pi)
plot(xlim,xlim + 2*pi)

subplot(2,1,2)
scatter(real(phsDiff),imag(phsDiff),[],spontFreq,'filled')
samelims
set(gca,'XLim',[-1, 1])
set(gca,'YLim',[-1, 1])
hold on
plot(real(phsors),imag(phsors))
grid on

c = colorbar;
c.Label.String = 'mean spontaneous frequency (Hz)';

% %% a plot
% figure; 
% cm = jet(length(unique(spontFreq)));
% for ii = 1:length(inboth)
%     arc = spontMag(ii).*exp(1i*linspace(spontPhs(ii),spontPhs(ii)+dtheta(ii),100));
%     plot(real(arc),imag(arc),'linewidth',2*spontMag(ii))
%     hold on
%     scatter(real(exp(1i*spontPhs(ii))),imag(exp(1i*spontPhs(ii))),'filled')
%     scatter(real(exp(1i*optoPhs(ii))),imag(exp(1i*optoPhs(ii))))
% end

%% per-event analysis 



