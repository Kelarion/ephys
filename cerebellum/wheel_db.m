n = 1;

db(n).name          = 'muad_dib';
db(n).serial        = 'M385799';
db(n).date          = '6-19-17';
db(n).depth         = {'2140','2217','2411','2579'};
db(n).side          = 'right';
db(n).badAux        = true;
db(n).hasOpto       = [true,true,true,true];
n = n + 1;

db(n).name          = 'leto_II';
db(n).serial        = 'M378800';
db(n).date          = '7-10-17';
db(n).depth         = {'2225','2525'};
db(n).side          = 'left';
db(n).badAux        = false;
db(n).hasOpto       = [true,true];
n = n + 1;

db(n).name          = 'ghanima';
db(n).serial        = 'M385800';
db(n).date          = '7-8-17';
db(n).depth         = {'2181','2278','2381'};
db(n).side          = 'right';
db(n).badAux        = false;
db(n).hasOpto       = [false,false,false];
n = n + 1;

db(n).name          = 'wensicia';
db(n).serial        = 'M38112';
db(n).date          = '8-28-17';
db(n).depth         = {'2331'};
db(n).side          = 'left';
db(n).badAux        = false;
db(n).hasOpto       = [true,true];
n = n + 1;

% db(n).name          = 'leto_II'; % the imec recordings
% db(n).serial        = 'M378800';
% db(n).date          = '7-11-17';
% db(n).depth         = {'3547'};
% db(n).side          = 'right';
% n = n + 1;