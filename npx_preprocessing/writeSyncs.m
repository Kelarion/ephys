% Script for aligning events and multiple probes to a common timeframe
% Taken from Sylvia, who I think got it from Nick, but made to work on
% batches of data.

%% all alignments (more at bottom)
if ~exist('overwrite','var'), overwrite = false; end
for k = 1:length(db)
    subject = db(k).mouse_name;
    date = db(k).date;
    
    %% Define folders
    subjectsFolder = '\\zubjects.cortexlab.net\Subjects';
    
    %% load the dataset
    [tags, hasEphys] = getEphysTags(subject, date);
    % sp = loadKSdir(subject, date);
    
    alignDir = fullfile(subjectsFolder, subject, date, 'alignments');
    if ~exist(alignDir, 'dir')
        mkdir(alignDir);
    elseif ~isempty(dir([alignDir '\correct_ephys*'])) && (~overwrite)
        continue
    end
    
    % determine what exp nums exist
    [expNums, blocks, hasBlock, pars, isMpep, tl, hasTimeline] = ...
        dat.whichExpNums(subject, date);
    TLexp = find(hasTimeline);
    
    %% plot drift map
    % incl = find(sp.spikeAmps > 40);
    % skip = round(length(incl)/1000000);
    % incl = incl(1:skip:end);
    % figure('Position', [1 33 1175 1083])
    % plotDriftmap(sp.st(incl), sp.spikeAmps(incl), sp.spikeDepths(incl));
    % axis tight
    % set(gca, 'box', 'off')
    
    %% align times (timeline to ephys)
    useFlipper = true;
    
    % for any ephys, load the sync data
    if hasEphys
        for t = 1:length(tags)
            if isempty(tags{t})
                [~, pdFlips, allET] = loadSync(subject, date);
            else
                [~, pdFlips, allET] = loadSync(subject, date, tags{t});
            end
            if useFlipper
                ephysFlips{t} = allET{5}{1};
            else
                ephysFlips{t} = pdFlips;
            end
        end
    end
    
    %% synchronize multiple ephys to each other
    if hasEphys
        if length(tags)>1
            for t2 = 2:length(tags)
                fprintf(1, 'correct ephys %s to %s\n', tags{t2}, tags{1});
                [~, b] = makeCorrection(ephysFlips{1}, ephysFlips{t2}, false);
                writeNPY(b, fullfile(alignDir, sprintf('correct_ephys_%s_to_ephys_%s.npy', tags{t2}, tags{1})));
            end
        end
    end
    
    %% detect sync events from timelines
    tlFlips = {};
    for e = 1:length(expNums)
        if hasTimeline(e)
            Timeline = tl{e};
            tt = Timeline.rawDAQTimestamps;
            if useFlipper
                evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'flipper'));
                evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
            else
                evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                evT = schmittTimes(tt, evTrace, [3 4]); % all flips, both up and down
                %         pdT = schmittTimes(tt, pd, [1.5 2]); % tried using TTL levels (0.8,2)
            end
            tlFlips{e} = evT;
        end
    end
    
    %% match up ephys and timeline events:
    % algorithm here is to go through each timeline available, figure out
    % whether the events in timeline align with any of those in the ephys. If
    % so, we have a conversion of events in that timeline into ephys
    %
    % Only align to the first ephys recording, since the other ones are aligned
    % to that
    if hasEphys
        ef = ephysFlips{1};
        if useFlipper && ef(1)<0.001
            % this happens when the flipper was in the high state to begin with
            % - a turning on event is registered at the first sample. But here
            % we need to drop it.
            ef = ef(2:end);
        end
        for e = 1:length(expNums)
            if hasTimeline(e)
                fprintf('trying to correct timeline %d to ephys\n', e);
                %Timeline = tl{e};
                tlT = tlFlips{e};
                
                success=false;
                if length(tlT)==length(ef)
                    % easy case: the two are exactly coextensive
                    [~,b] = makeCorrection(ef, tlT, false);
                    success = true;
                end
                if length(tlT)<length(ef) && ~isempty(tlT)
                    [~,b,success] = findCorrection(ef, tlT, false);
                end
                if success
                    writeNPY(b, fullfile(alignDir, ...
                        sprintf('correct_timeline_%d_to_ephys_%s.npy', ...
                        e, tags{1})));
                    fprintf('success\n');
                else
                    fprintf('could not correct timeline to ephys\n');
                end
            end
        end
    end
    
    %% match up blocks and mpeps to timeline in order:
    % want to connect each block or mpep with part of a timeline. So go through
    % each of these in order, looking through the timelines in sequence (of
    % what hasn't already been matched) looking for a match.
    lastTimes = zeros(1,length(expNums));
    for e = 1:length(expNums)
        if hasBlock(e)
            for eTL = 1:length(expNums)
                if hasTimeline(eTL)
                    fprintf('trying to correct block %d to timeline %d\n', e, eTL);
                    if useFlipper
                        % didn't get photodiode flips above, so get them now
                        Timeline = tl{eTL};
                        tt = Timeline.rawDAQTimestamps;
                        evTrace = Timeline.rawDAQData(:, strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                        pdT = schmittTimes(tt, evTrace, [3 4]);
                    else
                        pdT = tlFlips{eTL};
                    end
                    block = blocks{e};
                    sw = block.stimWindowUpdateTimes;
                    %                 sw = sw(2:end); % sometimes need this? Why? how did sw
                    % get an extra event at the beginning?
                    
                    success = false;
                    if length(sw)<=length(pdT) && length(sw)>1
                        [~,b,success,actualTimes] = findCorrection(pdT, sw, false);
                    end
                    if success
                        writeNPY(b, fullfile(alignDir, ...
                            sprintf('correct_block_%d_to_timeline_%d.npy', ...
                            e, eTL)));
                        writeNPY(actualTimes, fullfile(alignDir, ...
                            sprintf('block_%d_sw_in_timeline_%d.npy', ...
                            e, eTL)));
                        fprintf('  success\n');
                        lastTimes(eTL) = actualTimes(end);
                    else
                        fprintf('  could not correct block %d to timeline %d\n', e, eTL);
                    end
                end
            end
        elseif isMpep(e)
            for eTL = 1:length(expNums)
                if hasTimeline(eTL)
                    fprintf('trying to correct mpep %d to timeline %d\n', e, eTL);
                    p = pars{e}.Protocol;
                    nStims = numel(p.seqnums);
                    % An mpep stimulus has constant flips, with a gap in between
                    % stimuli. We'd like to know how many stimuli were shown, how long
                    % each lasted, and how much gap there was in between. But we can
                    % only really get the number of stimuli.
                    %                 minBetween = 0.2; % sec, assume at least this much time in between stimuli
                    %                 maxBetweenFlips = 2/60;
                    %                 if any(strcmp(p.parnames, 'dur'))
                    %                     estimatedDur = min(p.pars(strcmp(p.parnames, 'dur'),:))/10; % sec
                    %                     minDur = estimatedDur*0.75;
                    %                 else
                    %                     estimatedDur = [];
                    %                     minDur = 0.2;
                    %                 end
                    %                 nStims = numel(p.seqnums);
                    %                 pdT = tlFlips{eTL};
                    %                 pdT = pdT(pdT>lastTimes(eTL));
                    %
                    %                 dpdt = diff([0; pdT]);
                    %
                    %                 possibleStarts = find(dpdt>minBetween);
                    %                 possibleEnds = [possibleStarts(2:end)-1; length(pdT)];
                    %                 durs = pdT(possibleEnds)-pdT(possibleStarts);
                    %
                    %                 % gaps will be the gap *after* a stimulus
                    %                 gaps = pdT(possibleStarts(2:end))-pdT(possibleEnds(1:end-1));
                    
                    % dang, need a better system for mpep. The problem in
                    % Radnitz 2017-01-13 is that the photodiode was picking up
                    % some stuff around the sync square, which led to it being
                    % too bright in the down state to flip down for just one
                    % stimulus. Unfortunately that screws the whole thing.
                    Timeline = tl{eTL}; Fs = Timeline.hw.daqSampleRate;
                    tt = Timeline.rawDAQTimestamps;
                    tpd = Timeline.rawDAQData(:,strcmp({Timeline.hw.inputs.name}, 'photoDiode'));
                    
                    sig = conv(diff([0; tpd]).^2, ones(1,16/1000*Fs), 'same');
                    figure;
                    plot(tt(tt>lastTimes(eTL)), sig(tt>lastTimes(eTL)));
                    title(sprintf('expect %d stims', nStims));
                    
                    mpepStart = input('start of mpep? ');
                    mpepEnd = input('end of mpep? ');
                    thresh = input('lower/upper thresh? ');
                    
                    [flipTimes, flipsUp, flipsDown] = schmittTimes(tt, sig, thresh);
                    
                    flipsUp = flipsUp(flipsUp>mpepStart & flipsUp<mpepEnd);
                    flipsDown = flipsDown(flipsDown>mpepStart & flipsDown<mpepEnd);
                    
                    skippedFrames = (flipsUp(2:end)-flipsDown(1:end-1))<0.05; % assume stimuli are longer than 100ms
                    flipsUp = flipsUp(~[false; skippedFrames]);
                    flipsDown = flipsDown(~[skippedFrames; false]);
                    
                    if nStims==length(flipsUp)
                        fprintf(1, 'success\n');
                        stimOnsets = flipsUp;
                        stimOffsets = flipsDown;
                        
                        writeNPY(stimOnsets, fullfile(alignDir, ...
                            sprintf('mpep_%d_onsets_in_timeline_%d.npy', ...
                            e, eTL)));
                        writeNPY(stimOffsets, fullfile(alignDir, ...
                            sprintf('mpep_%d_offsets_in_timeline_%d.npy', ...
                            e, eTL)));
                        
                        lastTimes(eTL) = stimOffsets(end);
                    else
                        fprintf(1, 'unsuccessful - found %d, expected %d\n', length(flipsUp), nStims);
                    end
                end
            end
        end
    end
end
%% align eye video
% alignVideo(subject, date, TLexp, 'eye')
% resulting times are stored in Zserver\Data\EyeCamera; times need to be
% aligned to ephys times: eyeTimes.*timelineToEphys(1)+timelineToEphys(2)

%% Plot photodiode in ephys time to determine times of darkness
% bTLtoMaster = readNPY(fullfile(alignDir, ...
%     sprintf('correct_timeline_%d_to_ephys_%s.npy', TLexp, tags{1})));
% photodiode = tl{1}.rawDAQData(:,strcmp({tl{1}.hw.inputs.name}, 'photoDiode'));
% tlTime = tl{1}.rawDAQTimestamps;
% t = applyCorrection(tlTime, bTLtoMaster);
% figure('Position', [2 678 1917 420])
% plot(t, photodiode, 'k')
% xlabel('Time (aligned to ephys)')
% ylabel('photodiode')
% xlim([0 t(end)])
