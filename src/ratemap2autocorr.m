%{
    Create autocorrelation plot

    Author: Nate Sutton 2024
%}

%% Required Parameters %%
% specify rate map file to convert 
rm_filepath = "heat_maps_sim/heat_map_sml.mat";
% specify autocorrelogram file
ac_filepath = "saved_results/heat_map_sml_ac.mat";
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Autocorrelogram Generation %%
load(rm_filepath);
addpath cmbhome_utils
root = []; cell_selection=[1,1]; spk_x = []; spk_y = [];
heat_map_ac=plot_rate_map_ac(root, [1,1], heat_map, spk_x, spk_y);
colormap default
save(ac_filepath,"heat_map_ac");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%