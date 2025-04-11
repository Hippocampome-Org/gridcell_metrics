%%
%% Import CARLsim spikes using methods
%% from the CMBHome software
%%

function convert_carlsim_spikes(oat_location, x_positions_filepath, y_positions_filepath, spikes_binary_filepath, spikes_output_filepath, spk_bin_size, sel_nrn)
	%oat_location = "/comp_neuro/Software/CARLsim4_dgx_hc_09_18_21/tools/offline_analysis_toolbox"; % location of CARLsim's Offline Analysis Toolbox folder
	%x_positions_filepath = "/home/nmsutton/Downloads/temp/carlsim_results/results4/spikes/highres_pos_x.csv"; % file with virtual animal x-axis positions
	%y_positions_filepath = "/home/nmsutton/Downloads/temp/carlsim_results/results4/spikes/highres_pos_y.csv"; % file with virtual animal y-axis positions
	%spikes_binary_filepath = "/home/nmsutton/Downloads/temp/carlsim_results/results4/results/spk_MEC_LII_Stellate.dat"; % CARLsim's spikes data
	%spikes_output_filepath = "spikes.csv"; % conversion of spikes data into a csv format
	%spk_bin_size = 10; % spike reader bin size. Note: small bin sizes may take long processing with large spike sets. 40min sim with bin size 1 can take 10min to create plot.
	%sel_nrn=506; % neuron selected to analyze
	spk_bin_size = str2num(spk_bin_size);
	sel_nrn = str2num(sel_nrn);

	addpath(oat_location)
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