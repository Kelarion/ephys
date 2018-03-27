%% load 
load('C:\Users\Matteo\Documents\UNI\University of College\Ken class\allrips.mat')

% To convert from the 3-number identification (which correspond to cell ID,
% tag, and dataset) to a single unique cell ID, I use the below mappings. 
% They're equivalent to a conversion between base 5 and 10. For some 
% reason, using base-5 makes this mapping always invertible. (I checked 
% this with randomly-generated ID matrices where the 'digits' could vary 
% between 1 and 2000.)
wc2id = @(WC,b) WC*[b^2 b^1 b^0]';
id2wc = @(id,b) [floor(id/(b^2)) mod(floor(id/b),b) mod(id,b)];

CID = wc2id(R.whichCell,5);

