n = 0;

n = n+1;
db(n).mouse_name    = 'SS087';
db(n).date          = '2017-12-12';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2;
db(n).cwExp         = 3;
db(n).noiseExp      = 4;
db(n).passiveExp    = 5;
db(n).tags          = {'ZO' 'K2'};
% db(n).tagLocation   = {'' ''}; % either LA, LP, RA or RP 

n = n+1;
db(n).mouse_name    = 'SS087';
db(n).date          = '2017-12-13';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2; 
db(n).cwExp         = 3;
db(n).noiseExp      = 4;
db(n).passiveExp    = 5;
db(n).tags          = {'K1' 'K2' 'K3'};

n = n+1;
db(n).mouse_name    = 'SS087';
db(n).date          = '2017-12-14';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 3; 
db(n).cwExp         = 4;
db(n).noiseExp      = 5;
db(n).passiveExp    = 6;
db(n).tags          = {'K1' 'K2' 'K3'};

n = n+1;
db(n).mouse_name    = 'SS087';
db(n).date          = '2017-12-15';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2; 
db(n).cwExp         = 3;
db(n).noiseExp      = 4;
db(n).passiveExp    = 5;
db(n).tags          = {'K1' 'K2' 'K3'};

n = n+1;
db(n).mouse_name    = 'SS088';
db(n).date          = '2018-01-30';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2; 
db(n).cwExp         = 5;
db(n).noiseExp      = 6;
db(n).passiveExp    = 8;
db(n).tags          = {'K1','K2','K3','ZO'};

n = n+1;
db(n).mouse_name    = 'SS088';
db(n).date          = '2018-01-31';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2; 
db(n).cwExp         = 3;
db(n).noiseExp      = 4;
db(n).passiveExp    = 5;
db(n).tags          = {'K1','K2','K3','ZO'};

n = n+1;
db(n).mouse_name    = 'SS088';
db(n).date          = '2018-02-01';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2; 
db(n).cwExp         = 3;
db(n).noiseExp      = 4;
db(n).passiveExp    = 5;
db(n).tags          = {'K1','K2','ZO'};

n = n+1;
db(n).mouse_name    = 'SS088';
db(n).date          = '2018-02-02';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2; 
db(n).cwExp         = 3;
db(n).noiseExp      = 4;
db(n).passiveExp    = 5;
db(n).tags          = {'K1','K2','ZO'};

n = n+1;
db(n).mouse_name    = 'SS088';
db(n).date          = '2018-02-03';
db(n).dataServer    = 'main'; % means "zubjects\Subjects"
db(n).tlExp         = 2; 
db(n).cwExp         = 3;
db(n).noiseExp      = 4;
db(n).passiveExp    = 5;
db(n).tags          = {'K1','K2','K3','ZO'};


%% add nick's to the pile
r = findRecsWithArea('SCs');
q = findRecsWithArea('SCm');
p = dat.paths();
tmp = [];
for n = 1:length(r)
    j = strcmp({q.mouseName},r(n).mouseName) & strcmp({q.thisDate},r(n).thisDate);
    if any(j)
        db(end+1).mouse_name = q(j).mouseName;
        db(end).date = q(j).thisDate;
        db(end).dataServer = 'old';
        db(end).tlExp = q(j).tlExpNum;
        db(end).cwExp = q(j).cwExpNum;
        db(end).passiveExp = q(j).passiveExpNum;
%         db(end).noiseExp = q(j).noiseExpNum;
        db(end).tags = {q(j).tag};
    end
end