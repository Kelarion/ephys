% function extractSpkPara(input, kdir_ext, FileSave)
% input:
%   animal     = input.animal;
%   series     = input.iseries;
%   identifier = input.identifier;
%   probe      = input.probe;  (Default: Buzsaki32)
%   shanks     = input.shanks; (Default: shanks = 1:4)
%   PlotOn     = input.PlotOn; (Default: 0 i-> off)
%
% kdir_ext   = 1; klusta files saved in _klustakwik dir; 0 otherwise
% FileSave   = 1; saved output to *_info_*; 0 = off
%
% extract additional statistics from mean spk WF:
% minerr; FilAmp; PtT; PtP; FWHM; FW3M; Grad; Ratio; VarCh;
%
% Contributed by I-Chun Lin
% 2017-11-21 add two more metrics related to gradients (see lines marked %*)

function extractSpkPara(input, kdir_ext, FileSave)

if nargin < 3, FileSave = 1; end
if nargin < 2, kdir_ext = 1; end
if nargin < 1, warning('NOT ENOUGH INFORMATION'); end

code = '/mnt/zserver/Code/'; data = '/mnt/zserver/Data/';
addpath('/home/ichun/matlab/lsh');

animal     = input.animal;
series     = input.iseries;
identifier = input.identifier;
filename=[animal,'_s',num2str(series), '_', identifier];

kdir = [data,'multichanspikes/',animal,'/',num2str(series),'/'];
if kdir_ext>0, kdir = [kdir,'_klustakwik/']; end

if isfield(input, 'probe'),  probe  = input.probe; end
if isfield(input, 'shanks'), shanks = input.shanks; else shanks = 1:4; end

if ~isfield(input,'PlotOn'), PlotOn = 0; else PlotOn = input.PlotOn; end
if PlotOn>0
    load([data, 'multichanspikes/', animal, '/',num2str(series), '/', filename]);
    totChannel = length(SELECTED_CHANNELS);
    clear CHANNELS_ORDER inputStruct INTERNAL_REF_CONTACT lims ns5files2unify ...
        SELECTED_EXPERIMENTS SELECTED_CHANNELS;
end

for nsh = shanks % loop through shanks
    
    nElectrodes = 8; begChannel = nElectrodes*(nsh-1); clear excluCh;
    if exist('probe','var')
        if strcmp(probe,'edge32'), begChannel = 0; nElectrodes = 32; end
        if strcmp(probe,'edge16'), begChannel = 0; nElectrodes = 16; end
    end
    if isfield(input,'excluCh') & length(input.excluCh)>=nsh
        if length(input.excluCh)>0, excluCh = input.excluCh{nsh}; end
    end
    
    SpkTimes = double(hdf5read([kdir, filename, '.kwik'], ['channel_groups/', num2str(nsh-1), '/spikes/time_samples'])); nspikes = length(SpkTimes); % in sample #
    FilWFs = single(hdf5read([kdir,filename,'.kwx'], ['channel_groups/',num2str(nsh-1),'/waveforms_filtered']));
    load([filename, '_mSpk_sh',num2str(nsh-1)])
    load([filename, '_mFilSpk_sh',num2str(nsh-1)])
    load([filename, '_info_sh',num2str(nsh-1)])
    if PlotOn>0, plot_list = randperm(nspikes); plot_list = plot_list(1:60); end
    
    clear iCh nn_dist AvgWF_id Amp; % nn_idx
    iCh      = info.iCh;      % channel of the biggest amplitude
    AvgWF_id = info.AvgWF_id; % target & 4 NN same avg WF if good lsh dist
    nn_dist  = info.nn_dist;  % 1st, 4th & every 10th nn
    Amp      = info.Amp';     % amplitude from avg detrened WF
    
    clear lll ll nlsh index; ll = int8(zeros(nspikes, 1));
    lll = find((double(AvgWF_id)-[1:nspikes])~=0); ll(lll) = 1;
    for i=1:nspikes, nlsh(i,1)=14-length(find(isnan(nn_dist(:,i)))); end
    index = find(info.AvgWF_id'~=0 & info.iCh>0 & (nlsh~=0 | ll~=0));
    
    minerr = single(nan(nspikes, 1));   % fitting error
    FilAmp = single(nan(nspikes, 1));   % amplitude from filtered avg detrended WF
    Grad   = single(nan(nspikes, 2));   % gradient @.067+0.5 ms
    VarCh  = single(nan(nspikes, 1));   % variablity across channels
    PtT    = int8(zeros(nspikes, 1));   % trough position
    PtP    = int8(zeros(nspikes, 1));   % peak position
    Ratio  = single(zeros(nspikes, 1)); % Height ratio (peak/trough)
    FWHM   = single(zeros(nspikes, 1)); % full width 1/2 maximum
    FW3M   = single(zeros(nspikes, 1)); % full width 1/3 maximum
    PreMax = single(zeros(nspikes, 1));
    
    for ii = 1:length(index), ispk = index(ii);
        clear aa dTrdAvgFilWF dTrdAvgWF maxch amp tmp1 tmp2 t1 y1 t y c d H1 H2;
        if rem(ispk, 1e4)==0, txt = sprintf('ispk = %d', ispk); disp(txt); end
        dTrdAvgFilWF = dTrdAvgFilWFs(:,:,AvgWF_id(ispk)); % detrended avg filtered WF
        dTrdAvgWF    = dTrdAvgWFs(:,:,AvgWF_id(ispk));    % detrended avg raw WF
        maxch        = iCh(ispk); % channel of max amplitude
        amp          = Amp(ispk); % amplitude (trough)
        
        % peak and trough positions, gradients @ 0.067 & 0.5 ms from trough, FWXM
        [vt, c]  = min(dTrdAvgWF(maxch,17:23 )); PtT(ispk) = 16 + c;
        tmp1 = smooth(dTrdAvgWF(maxch,:)', 0.1, 'loess');
        tmp2 = smooth(dTrdAvgWF(maxch,:)', 0.5, 'loess');
        [vp, c]  = max(tmp1(PtT(ispk)+1:45)); PtP(ispk) = PtT(ispk) + c;
        [vp1, PtP1] = max(tmp2(PtT(ispk)+1:45)); PtP1 = PtT(ispk) + PtP1;
        if PtP1<PtP(ispk), PtP(ispk) = PtP1; vp = vp1; end
        aa = diff(tmp1); PtP1 = min(cat(1,99,find(aa(PtT(ispk)+3:end)<0)))+PtT(ispk);
        if PtP1<PtP(ispk), PtP(ispk) = PtP1; vp = tmp1(PtP1); end
        Grad(ispk,1) = aa(PtT(ispk)+2); Grad(ispk,2) = aa(PtT(ispk)+15);
        Ratio(ispk)  = -vp/vt;
        
        t1 = single(PtT(ispk)-5:PtT(ispk)+9); y1 = tmp1(t1)'; t = t1(1):.1:t1(end);
        y = interp1(single(t1), y1, t); d = find(t==PtT(ispk));
        [~, H1] = min(abs(y(1:d)+amp/2));   [~, H2] = min(abs(y(d:end)+amp/2));
        H1 = t(H1); H2 = t(H2+d-1); FWHM(ispk) = H2-H1;
        [~, H1] = min(abs(y(1:d)+2*amp/3)); [~, H2] = min(abs(y(d:end)+2*amp/3));
        H1 = t(H1); H2 = t(H2+d-1); FW3M(ispk) = H2-H1;
        
        % variabilty across channels, fitting error, filtered amplitude
        if exist('excluCh','var'), clear ble;
            ble = dTrdAvgFilWF; ble(excluCh,:,:) = []; VarCh(ispk) = max(std(ble));
        else VarCh(ispk) = max(std(dTrdAvgFilWF));
        end
        dTrdAvgFilWF = dTrdAvgFilWF(maxch,:);
        [err curves] = min_err(FilWFs(maxch,:,ispk)', dTrdAvgFilWF');
        [minerr(ispk), bestPt] = min(err); curve1 = curves(bestPt,:);
        FilAmp(ispk)  = max(dTrdAvgFilWF(1:10))-min(dTrdAvgFilWF(7:13));
        
        if PlotOn>0 & ismember(ii, plot_list), clear TmpWFs;
            figure; subplot(3,3,1:2);
            plot([1:4 10:10:100],nn_dist(:,AvgWF_id(ispk)),'o-'); xlim([0 100]); grid on
            title(gca, sprintf('ID=%d (%d) amp famp=%.0f %.0f var=%.0f peak=%d grad1 grad2=%.1f %.1f ratio=%.2f', ispk, AvgWF_id(ispk), amp, FilAmp(ispk), VarCh(ispk), PtP(ispk), Grad(ispk,1), Grad(ispk,2), Ratio(ispk)))
            subplot(3,3,3); plot(dTrdAvgWF(maxch,:)', 'ko-'); hold on; grid on;
            plot(tmp1, 'r-'); plot(tmp2, 'm-'); plot(aa, 'x-'); clear aa;
            subplot(3,3,4:6);
            plot([1:4 10:10:100], nn_dist(:,AvgWF_id(ispk))/FilAmp(ispk),'o-'); hold on
            plot([0 100], [0.025 0.025], 'r--'); xlim([1 100]); grid on
            
            TmpWFs = single(SpikeExtract([filename,'.dat'], SpkTimes(ispk), totChannel, begChannel+1:begChannel+nElectrodes, 20, 50));
            subplot(3,3,7); plot(TmpWFs'); hold on;
            plot(TmpWFs(maxch,:), 'o'); xlim([1 70]); grid on;
            subplot(3,3,8); plot(dTrdAvgWF'); hold on;
            plot(dTrdAvgWF(maxch,:)','o'); xlim([1 70]); grid on;
            subplot(3,3,9); plot(dTrdAvgFilWF', 'ko-'); hold on; grid on; xlim([1 70]);
            aa=1+(4-bestPt):20+(4-bestPt); beg_aa=find(aa==1); aa=aa(find(aa>0 & aa<20));
            if length(beg_aa)<1 | beg_aa<1, plot(curve1(aa), 'r-');
            else plot(beg_aa:beg_aa+length(aa)-1, curve1(aa), 'r-');
            end
            title(gca,sprintf('%.2f %.1f %.1f', minerr(ispk)/FilAmp(ispk), FWHM(ispk), FW3M(ispk)));
            savefig(['spk', num2str(ispk), '_', filename, '.fig']); close all;
        end % if PlotOn
        
    end % for ii
    
    info.minerr = minerr; info.FilAmp = FilAmp;
    info.PtT    = PtT;    info.PtP    = PtP;
    info.FWHM   = FWHM;   info.FW3M   = FW3M;
    info.Grad   = Grad;   info.Ratio  = Ratio;
    info.VarCh  = VarCh;  info.PreMax = PreMax;
    
    if length(find(PtP~=0))==length(index) & FileSave > 0
        save([filename, '_info_sh', num2str(nsh-1)], 'info');
    else warning('Post-analysis did not complete, temp info file saved');
        save(['Tmp_', filename, '_info_sh', num2str(nsh-1)], 'info');
    end
    
end % for nsh = shanks

