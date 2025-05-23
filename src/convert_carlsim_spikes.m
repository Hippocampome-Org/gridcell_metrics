%%
%% Import CARLsim spikes using methods
%% from the CMBHome software
%%

function convert_carlsim_spikes(oat_location, x_positions_filepath, y_positions_filepath, spikes_binary_filepath, spikes_output_filepath, spk_bin_size, sel_nrn)
	spk_bin_size = str2num(spk_bin_size);
	sel_nrn = str2num(sel_nrn);

    x = readmatrix(x_positions_filepath);
    y = readmatrix(y_positions_filepath);
	spk_data = SpikeReader(spikes_binary_filepath, false, 'silent');
	spikes = spk_data.readSpikes(spk_bin_size);
	spikes_output_file = fopen(spikes_output_filepath,'w');

    spk_t=find(spikes(:,sel_nrn)~=0);
    spk_t=spk_t*spk_bin_size;
    spk_x = []; spk_y = [];

    for i=1:length(spk_t)
        spk_x=[spk_x,x(spk_t(i))];
        spk_y=[spk_y,y(spk_t(i))];
    	fprintf(spikes_output_file,"%f,%f,%f\n",spk_t(i),spk_x(i),spk_y(i));
    end

    fclose(spikes_output_file);
end