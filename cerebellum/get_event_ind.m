function [ i_cross ] = get_event_ind( fname,nchan,chan2scan,thresh,invert )
% Get indices of digital pulses from binary analog data file. If invert = 1, it
% looks for negative crossings.
% [ i_cross ] = get_event_ind( fname,nchan,chan2scan,thresh,invert )


i_step = 50000*60*5;
%scale = 2.5./(200*(2^15));
f = fopen(fname);
%olap = 50000;
%thresh = .001;
if isempty(thresh)
    thresh = nan(size(chan2scan));
end

% Iterate over blocks.
fseek(f,0,'eof');
sz = ftell(f)/(2*nchan);
frewind(f);
n_block = ceil(sz/i_step);
i_cross = cell(1,length(chan2scan));

rawdat = memmapfile(fname, 'Format',{'int16', [nchan sz], 'x'});
for i = 1:n_block
    
    disp(['Processing block ' num2str(i) '/' num2str(n_block) '... ']);
    clear data_tmp
    
    % Load data.
    i_offset = (i-1)*i_step;
    
%     for j = 1:length(chan2scan)
%         data_tmp(j,:) = load_channel('int16',fname,nchan,chan2scan(j),i_offset,i_step);
%         if invert == 1
%             data_tmp(j,:) = -data_tmp(j,:);
%         end
%     end
    
    % Find pulses.
    for j = 1:length(chan2scan)
        blockSamps = i_offset + [1:i_step];
        data_tmp = rawdat.Data.x(chan2scan(j),blockSamps(blockSamps<sz));
        %{
        [amp_tmp i_cross_tmp] = findpeaks(diff(data_tmp(j,:)),'minpeakheight',thresh(j));
        if ~isempty(i_cross_tmp)
            i_cross{j} = [i_cross{j} i_cross_tmp+(i-1)*i_step];
            amp{j} = amp_tmp;
        end
        %}
        if isnan(thresh(j)) % automatic threshold setting (Matteo)
            gm = fitgmdist(double(data_tmp(1:100:end))',2);
            %                         x = linspace(gm.mu(1),gm.mu(2),100);
            %                         sig = pdf('norm',x,gm.mu(1),gm.Sigma(1));
            %                         noys = pdf('norm',x,gm.mu(2),gm.Sigma(2));
            %             [~,cross] = min(abs(sig./noys-1));
            %             thresh(j) = abs(x(cross) - gm.mu);
            %             thresh(j) = abs(diff(gm.mu)/2);
            thresh(j) = mean(gm.mu);
        end
        
        %         i_cross_tmp = find(abs(diff(data_tmp(j,:))) > thresh(j));
        %         i_cross_tmp(diff([0 i_cross_tmp])<100) = []; % different detection algo (Matteo)
        i_cross_tmp = find(data_tmp(1:(end-1))<thresh(j) & data_tmp(2:end)>=thresh(j));
        i_cross{j} = [i_cross{j} i_cross_tmp+i_offset];

    end
    
    
end

for i = 1:length(i_cross)
    i_cross{i} = sort(i_cross{i});
end

end