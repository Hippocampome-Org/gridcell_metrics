function smooth_rate_map(spikes_file, output_file, Xs_file, Ys_file, occupancy_norm, omit_islands, omit_noocc, binside, std_smooth_kernel, plot_subsect, grid_size, plot_size, plot_spikes, fs_video, timestep, use_smoothing, show_plot, save_file)
% spikes = readmatrix("input_data/spikes.csv");
% output_file = "saved_results/rate_map.txt";
% x = readmatrix("/media/nmsutton/StorageDrive/comp_neuro/gmu/research/ach_sim/data/neuron_23/highres_pos_x.csv");
% y = readmatrix("/media/nmsutton/StorageDrive/comp_neuro/gmu/research/ach_sim/data/neuron_23/highres_pos_y.csv");
% occupancy_norm = 1; % perform occupancy normalization. this adjusts rates by number of visits to locations.
% omit_islands = 1; % see later comments
% omit_noocc = 1; % set no occupancy to zero
% fs_video = 50; % sampling rate from video (samples/sec)
% binside = 3;
% std_smooth_kernel = 3.333;
% plot_subsect=1;
% grid_size=40;
% plot_size=31;
% plot_spikes = 1;
% timestep=20; % recordings timestep of samples in ms
% use_smoothing=1;
% show_plot=0;
% save_file=1;
spikes = readmatrix(spikes_file);
x = readmatrix(Xs_file);
y = readmatrix(Ys_file);
occupancy_norm = str2num(occupancy_norm);
omit_islands = str2num(omit_islands);
omit_noocc = str2num(omit_noocc);
binside = str2num(binside);
std_smooth_kernel = str2num(std_smooth_kernel);
plot_subsect = str2num(plot_subsect);
grid_size = str2num(grid_size);
plot_size = str2num(plot_size);
plot_spikes = str2num(plot_spikes);
fs_video = str2num(fs_video);
timestep = str2num(timestep);
use_smoothing = str2num(use_smoothing);
show_plot = str2num(show_plot);
save_file = str2num(save_file);

spikes_x = spikes(1:end,3);
spikes_y = spikes(1:end,2);

% update grid size if larger x or y values are found
if max(spikes_x)>grid_size || max(spikes_y)>grid_size
	if max(spikes_x)>max(spikes_y)
		grid_size=ceil(max(spikes_x));
        if floor(min(spikes_x)) < 0
            grid_size=ceil(max(spikes_x))-floor(min(spikes_x));
        end
	else
		grid_size=ceil(max(spikes_y));
        if floor(min(spikes_y)) < 0
            grid_size=ceil(max(spikes_y))-floor(min(spikes_y));
        end
	end
end
% shift matrix to account for negative values
if floor(min(spikes_x)) < 0
    spikes_x=spikes_x-floor(min(spikes_x));
end
if floor(min(spikes_y)) < 0
    spikes_y=spikes_y-floor(min(spikes_y));
end

heat_map = zeros(1,grid_size*grid_size);
xdim = linspace(0,grid_size-1,grid_size);
ydim = linspace(0,grid_size-1,grid_size);
if occupancy_norm
	occupancy = hist3([x,y],'Edges',{xdim, ydim})/fs_video; 
	no_occupancy = occupancy==0; % mark indeces where there was no occupancy so we can correct after smoothing
end
heat_map = hist3([spikes_x,spikes_y],'Edges',{xdim, ydim});
heat_map = heat_map*timestep;
if plot_subsect
    s = ceil(((grid_size-plot_size)/2));
    e = ceil(plot_size+s);
    if occupancy_norm
        occupancy = occupancy(s:e, s:e); % crop to intended plot size
        no_occupancy = no_occupancy(s:e, s:e); % crop to intended plot size
    end
    heat_map = heat_map(s:e, s:e); % crop to intended plot size
end
if omit_islands
	% Takes matrix of occupancy values, calculates center of mass of pixels>0,
	% and then finds all disconnected pixels and sets them = 0
	s = regionprops(occupancy>0, {'FilledArea', 'PixelIdxList'});
	l = numel(s);
	areas = vertcat(s.FilledArea);
	[~, inds] = sort(areas);
	for i = 1:length(inds)-1				    
	    occupancy(s(inds(i)).PixelIdxList) = 0;				    
	end
end
if occupancy_norm
	occupancy=occupancy';
end
if use_smoothing
	if occupancy_norm
		heat_map = SmoothMat(heat_map, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside)./SmoothMat(occupancy, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside);
	else
		heat_map = SmoothMat(heat_map, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside); % smooth the spikes and occupancy with a 5x5 bin gaussian with std=1
	end
end
if omit_noocc
    heat_map(no_occupancy) = 0; % set no occupancy to zero
end

if show_plot
    imagesc(heat_map);
    axis('tight')
    xlabel('animal position on x axis') 
    ylabel('animal position on y axis')
    cb = colorbar;
    ax = gca;
    ax.YDir = 'normal'; % have y-axis 0 on bottom left
end

if save_file
    writematrix(heat_map,output_file);
end
end