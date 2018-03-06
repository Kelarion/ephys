function [ data ] = load_channel( int_type,fname,nchan,chan2get,i_start,n_samp )
% [ data ] = load_bin_single_channel( int_type,fname,nchan,chan2get,i_start,n_samp )

%{
gain = 200;
bitdepth = str2num(int_type((end-1):end));
scale = 2.5./(gain*(2^(bitdepth-1)));
%}

f = fopen(fname);
byte_offset = i_start*nchan*2 + (chan2get-1)*2;
fseek(f,byte_offset,'bof');
data = fread(f,[1 n_samp],int_type,2*(nchan-1));
%data = scale.*data;
fclose(f);

end