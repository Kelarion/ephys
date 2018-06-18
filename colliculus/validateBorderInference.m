%% make db
localRoot = 'C:\DATA\Spikes\';

clear db
db = struct;
r = findRecsWithArea('SCs');
q = findRecsWithArea('SCm');
p = dat.paths();
ii = 1;
for n = 1:length(r)
    j = strcmp({q.mouseName},r(n).mouseName) & strcmp({q.thisDate},r(n).thisDate);
    if any(j)
        db(ii).mouse_name = q(j).mouseName;
        db(ii).date = q(j).thisDate;
%         db(ii).dataServer = q(j).whichServer;
        db(ii).tlExp = q(j).tlExpNum;
        db(ii).cwExp = q(j).cwExpNum;
        db(ii).passiveExp = q(j).passiveExpNum;
        db(ii).noiseExp = q(j).noiseExpNum;
        db(ii).ksRoot = [p.(['oldRepository']) '\'];
        db(ii).tags = {q(j).tag};
        ii = ii+1;
    end
end

visrecs = zeros(length(db),1);
q = findRecsWithArea('VISp');
for n = 1:length(q)
    j = strcmp({db.mouse_name},q(n).mouseName) & strcmp({db.date},q(n).thisDate);
    visrecs = visrecs | j(:);
end

%%
inftop = nan(length(db),1);
infbot = nan(length(db),1);
knowntop = nan(length(db),1);
knownbot = nan(length(db),1);
knownmid = nan(length(db),1);
for k = 1:length(db)
    for t = 1:length(db(k).tags), thisTag = db(k).tags{t};
        % load everything for that probe
        dsetFolders = [db(k).mouse_name '\' db(k).date '\'];
        
        [~,dataDir, ~, alfDir] = ... 
            expDirs(db(k).mouse_name,db(k).date,thisTag,'old');
        ksDir = [dataDir '\sorting\'];
        spks = loadNeuralData(ksDir);
        
        saveFolder = [localRoot dsetFolders 'ephys_' thisTag '\'];
        bfname =[alfDir thisTag '\borders_' thisTag '.tsv'];
        snname = [saveFolder 'sparse_noise_RFs.mat'];
        if exist(bfname,'file') && exist(snname,'file')
            disp('Reading histology data...'); 
            bord = readtable(bfname ,'Delimiter','\t', 'FileType', 'text');
            sc = contains(lower(bord.name),'superior colliculus');
            scupperGT{k} = bord.upperBorder(sc);
            sclowerGT{k} = bord.lowerBorder(sc);
            nyms{k} = bord.acronym(sc);
            
            [scschan, scmchan] = wheresCortex([dataDir spks.lfp_path],'snrf',snname);
            if isempty(scschan)
                t0 = aln.noise2tl(2) + aln.tl2ref(2) - aln.tag2ref(2);
                scschan = wheresCortex([dataDir spks.lfp_path],'corrs',[t0 t0+400]);
                scmchan = 0;
            end
            scTopInf{k} = spks.ycoords(scschan);
            scsBotInf{k} = spks.ycoords(scmchan);
            % package
            inftop(k) = scTopInf{k};
            infbot(k) = scsBotInf{k};
            knowntop(k) = max(scupperGT{k});
            knownbot(k) = min(sclowerGT{k});
            if any(contains(nyms{k},'sg'))
                knownmid(k) = sclowerGT{k}(contains(nyms{k},'sg'));
            end
            
        else
            continue
        end
        
    end
end

%% plot
figure('units','normalized','position',[0.1234 0.4241 0.8333 0.3852]);
subplot(1,3,1)
scatter(knowntop,inftop,'filled')
text(knowntop(visrecs)+10,inftop(visrecs)+10,'V1','fontsize',11)
samelims; axis square;
hold on; plot(xlim,xlim,'--k')
ylabel('inferred border (um)')
xlabel('known border (um)')
title('SC surface')

subplot(1,3,2)
scatter(knownmid,infbot,'filled')
text(knownmid(visrecs)+10,infbot(visrecs)+10,'V1','fontsize',11)
samelims; axis square;
hold on; plot(xlim,xlim,'--k')
xlabel('known border (um)')
title('SCs/SCm')

subplot(1,3,3)
scatter(knownbot,infbot,'filled')
text(knownbot(visrecs)+10,infbot(visrecs)+10,'V1','fontsize',11)
samelims; axis square;
hold on; plot(xlim,xlim,'--k')
xlabel('known border (um)')
title('SCm/other midbrain')

