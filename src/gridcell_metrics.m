%{
    Automatic grid cell metric reporter
    
    Author: Nate Sutton, 2023
    
    Required parameter:
    - nonfld_filt_perc: percentage of lowest firing rates to remove from
    data. This removes non-field firing.

    Optional parameters:
    - use_ac: choose to use autocorrelogram (ac) (1) or standard rate map (0)
    - only_center_seven: filter out fields except for seven closest to the plot center
    - dist_thresh: maximum field distance to be considered in group of
    field distances.
    - fld_sz_thresh: minimum field size that has contact with a border for
    the field to avoid being filtered out.
    - manual_field_exclud: field indices to exclude from spacing
    statistics.
    - load_px2cm_conv: load saved list of pixel to centimeter distance
    conversions.
    - spac_exclud; size_exclud; ang_exclud: specify any fields to exclude 
    from size, spacing, or angle statistics.
    dbscan_epsilon; dbscan_min_pts: db_scan parameters

    Plotting:
    - plot_fields_detected: plot the detected fields.
    - plot_orig_firing: plot the rate map or autocorrelogram used as neural
    data input.
    
    references: https://www.mathworks.com/help/stats/dbscan.html
    https://www.mathworks.com/matlabcentral/answers/501948-how-to-automatically-find-the-angle-from-2-points
    https://www.mathworks.com/matlabcentral/answers/629408-how-to-use-convolution-on-a-2d-matrix-with-a-kernel
    http://agamtyagi.blogspot.com/2013/02/matlab-code-for-famous-mexican-hat.html
    https://gamedev.stackexchange.com/questions/4467/comparing-angles-and-working-out-the-difference
%}

%% Parameter Options %%
%% Required Settings %%
px2cm=100/32; % conversion of each pixel to a centimeter lenth in an enviornment
nonfld_filt_perc=.31; % non-filtered percentage threshold
load_plot_from_file=1; % load neural data plot from individual file (1; plot_filepath) or instead pick a file from a numbered list (0; heat_maps_list or heat_maps_ac_list)
plot_filepath="heat_maps_real_ac/n20_ac.mat"; % path to file that contains the neural data
if load_plot_from_file==0
    load("heat_maps_list.mat"); small=35; medium=36; large=37; % list of rate map plots. real cell file_number 15 = small, 23 = medium, 17 = large
    load("heat_maps_ac_list.mat"); % list of autocorrelogram plots
    heat_map_selection=89; % pick plot number from the supplied list
end
%% Optional Settings %%
use_ac=1; % choose to use autocorrelogram (ac) (1) or standard rate map (0)
use_dist_thresh=0;
dist_thresh=500;
use_fld_sz_thresh=0;
fld_sz_thresh=6;
spac_exclud = []; size_exclud = []; ang_exclud = []; % specify any fields to exclude from statistics
manual_field_exclud = 1; % field indices to exclude from spacing statistics
if manual_field_exclud==1 excluded_fields; end
load_px2cm_conv = 1; % load saved list of pixel to centimeter distance conversions
if load_px2cm_conv==1 saved_px2cm_conv; end
dbscan_epsilon=3; dbscan_min_pts=1; % db_scan parameters
load_custom_rm=0; % load custom rate map saved from a file
load_custom_ac=0; % load custom ac saved from a file
only_center_seven=1; % filter out fields except for seven closest to the plot center
use_tophat_filter=0;
use_centsurr_filter=0;
if exist('heat_map_selection_custom','var') heat_map_selection=heat_map_selection_custom; end
if use_ac == 0 only_center_seven=0; end
minimal_plotting_mode=0;
if exist('minimal_plotting_mode_custom','var') minimal_plotting_mode=minimal_plotting_mode_custom; end
auto_export_plots=0;
if exist('auto_export_plots_custom','var') auto_export_plots=auto_export_plots_custom; end
filename_sizes='saved_results/field_size_records.txt';
if exist('filename_sizes_custom','var') filename_sizes=filename_sizes_custom; end
filename_spacings='saved_results/field_spacing_records.txt';
if exist('filename_spacings_custom','var') filename_spacings=filename_spacings_custom; end
filename_rotations='saved_results/save_field_rotation.txt';
if exist('filename_rotations_custom','var') filename_rotations=filename_rotations_custom; end
if exist('saved_nonfld_filt_perc','var') nonfld_filt_perc=saved_nonfld_filt_perc; end
% ac_xaxis_dim=63; % commented out because this may no longer be used. for custom ac, use this x-axis length. This is the length in pixels of the x- and y-axes of the autocorrelogram plot
%% Plotting Settings %%
plot_fields_detected=1;
plot_orig_firing=1; % plots the original firing data's rate map or autocorrelogram
if exist('plot_fields_detected_custom','var') plot_fields_detected=plot_fields_detected_custom; end
if exist('plot_orig_firing_custom','var') plot_orig_firing=plot_orig_firing_custom; end
plot_legend=0; print_angles=0;
control_window_size=0;
custom=0; custom2=0; custom3=0; custom4=0;
if custom load_custom_rm=1; plot_orig_firing=0; end
if custom2 load_custom_ac=1; plot_orig_firing=0; end
if custom3 load_custom_rm=1; use_ac=0; only_center_seven=0; end
if custom4 load_custom_rm=0; use_ac=1; only_center_seven=1; end
sml_ang_cnt_fld = 361; % smallest angle from center field. Default is 361 which is beyond 360 which could be a highest angle possible.
sml_ang_cent_num = 0; % centriod number of sml_ang_cnt_fld match
%%%%%%%%%%%%%%%%%%%%%%

%% Load Neural Data %%
if plot_fields_detected && plot_orig_firing control_window_size=1; end
if exist('use_newly_generated_rm','var')
    % only load heat map from file if no heat map is already in memory that
    % is indicated by use_newly_generated_rm
    if str2num([string(use_newly_generated_rm)]) == 1
        use_ac=0;
        only_center_seven=0;
    else
        load(heat_maps(heat_map_selection));
    end
elseif load_plot_from_file == 0
    load(heat_maps(heat_map_selection));
end
if load_custom_rm
    load("heat_maps_real_custom_list.mat");
    load(heat_maps_real_custom_list(heat_map_selection));
    heat_map=heat_map_custom;
end
if load_plot_from_file == 1
    load(plot_filepath);
    if exist('heat_map_ac','var') heat_map=heat_map_ac; end
end
if use_ac && load_plot_from_file==0
    load(heat_maps_ac_list(heat_map_selection));
    heat_map=heat_map_ac;
end
if load_custom_ac
    load("heat_maps_ac_custom_list.mat");
    load(heat_maps_ac_custom_list(heat_map_selection));
    heat_map=heat_map_custom;
end
if exist('use_newly_generated_ac','var')
    if str2num([string(create_ac)]) == 1
        heat_map=auto_corr_rm;
    end
end
heat_map_orig=heat_map;
res = size(heat_map,2); % x-axis resolution
heat_map2=[];
if load_custom_ac==0
    % rescale
    peak = max(max(heat_map));
    heat_map = heat_map * 1/peak;
    
    % filter non-fields
    heat_map(find(heat_map<nonfld_filt_perc))=0.0;

    % create filtered heat map for saving to a file
    heat_map_custom=heat_map;
    heat_map_custom(find(heat_map_custom>0))=1.1111;
    heat_map=heat_map_custom;
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Neural Data %%

if use_tophat_filter
    se = strel('disk',5);
    tophatFiltered = imtophat(heat_map_orig,se);
    heat_map=tophatFiltered;
end
if use_centsurr_filter
    centsurr_height=1;
    centsurr_width=1.05;
    % convolution with center-surround matrix
    % create center-surround matrix
    [x,y] = meshgrid(-8:0.5:8);
    %[x,y] = meshgrid(-15.5:0.5:15.5);
    r = sqrt(x.^2+y.^2)+eps;
    centsurr = (sin(r*centsurr_width)*centsurr_height)./r;
    heat_map = basic_convolution(heat_map_orig,centsurr);
    % plot center-surround matrix
    % surfc(centsurr), shading flat
    % axis off
end

% find fields
for y=1:res
    for x=1:res
        if heat_map(y,x) > 0.0
            heat_map2=[heat_map2;[y,x]];
        end
    end
end

idx = dbscan(heat_map2,dbscan_epsilon,dbscan_min_pts);

% sort fields
% organize field points into matrix rows
fields_num=max(idx);
points_per_field = histcounts(idx);
most_points=max(points_per_field);
fields_x=zeros(fields_num,most_points);
fields_y=zeros(fields_num,most_points);
for i=1:fields_num
    pnt_ctr=0;
    for j=1:size(idx,1)
        if idx(j)==i
            pnt_ctr=pnt_ctr+1;
            fields_x(i,pnt_ctr)=heat_map2(j,2);
            fields_y(i,pnt_ctr)=heat_map2(j,1);
        end
    end
end

% filter out fields at edges without enough points
filter_out=[];
if use_fld_sz_thresh==1
    for i=1:fields_num
        fprintf("processing: %d\n",i);
        non_zero=find(fields_x(i,:)~=0);
        non_zero_sz=size(non_zero,2);
        fld_x=[];fld_y=[];
        fld_x=fields_x(i,non_zero);
        fld_y=fields_y(i,non_zero);
        % search for border contact
        border_w=find(fld_x<2);
        border_e=find(fld_x>(res-1));
        border_s=find(fld_y<2);
        border_n=find(fld_y>(res-1));
        if isempty(border_n) border_n=0; else border_n=1; end
        if isempty(border_s) border_s=0; else border_s=1; end
        if isempty(border_e) border_e=0; else border_e=1; end
        if isempty(border_w) border_w=0; else border_w=1; end
        if border_w || border_n || border_e || border_s
            if non_zero_sz<=fld_sz_thresh
                filter_out=[filter_out,i];
            end
        end
    end
end

keep_in=[];
for i=1:fields_num
    if isempty(find(filter_out==i))
        keep_in=[keep_in,i];
    end
end
fields_num=fields_num-size(filter_out,2);
fields_x2=fields_x(keep_in,:);
fields_y2=fields_y(keep_in,:);
fields_x=fields_x2;
fields_y=fields_y2;

% find centroids
centroid_x=[];
centroid_y=[];
for i=1:fields_num
    sum_x=0;
    sum_y=0;
    pnt_ctr=0;
    for j=1:most_points
        if fields_x(i,j)~=0
            sum_x=sum_x+fields_x(i,j);
            sum_y=sum_y+fields_y(i,j);
            pnt_ctr=pnt_ctr+1;
        end
    end
    centroid_x=[centroid_x,sum_x/pnt_ctr];
    centroid_y=[centroid_y,sum_y/pnt_ctr];
end

% optional filter to only keep 7 fields closest to the center
closest_seven=[];
field_distances2=[];
if only_center_seven
    center_point=size(heat_map_ac,2)/2;
    for i=1:fields_num
        field_dist=euc_d(centroid_x(i),centroid_y(i),center_point,center_point);
        field_distances2=[field_distances2;[i,field_dist]];
    end
    closest_seven=sortrows(field_distances2,2);
    fld_cnt = size(field_distances2,1);
    if fld_cnt > 7
        fld_cnt = 7;
    end
    closest_seven=closest_seven(1:fld_cnt);
    % closest_seven=closest_seven(1:7);
    % filter out unwanted fields
    % filter_out=(fields_num-7);
    filter_out=(fields_num-fld_cnt);
    % fields_num=7;
    fields_num=fld_cnt;
    fields_x=fields_x(closest_seven,:);
    fields_y=fields_y(closest_seven,:);
    centroid_x=centroid_x(:,closest_seven);
    centroid_y=centroid_y(:,closest_seven);
    % only keep points in the seven fields
%     heat_map2=[];
%     idx=[];
%     for i=1:fields_num
%         for j=1:size(fields_y,2)
%             if fields_y(i,j)~= 0
%                 heat_map2=[heat_map2;[fields_y(i,j),fields_x(i,j)]];
%                 idx=[idx;i];
%             end
%         end
%     end
end

heat_map2=[];
idx=[];
for i=1:fields_num
    for j=1:size(fields_y,2)
        if fields_y(i,j)~= 0
            heat_map2=[heat_map2;[fields_y(i,j),fields_x(i,j)]];
            idx=[idx;i];
        end
    end
end

% find field spacing
field_distances=[];
% find most center field
center_field_idx=1; % center field index
center_x=size(heat_map,1)/2;
center_y=size(heat_map,2)/2;
smallest_dist=euc_d(centroid_x(1),centroid_y(1),center_x,center_y);
for i=1:fields_num
    field_dist=euc_d(centroid_x(i),centroid_y(i),center_x,center_y);
    if (field_dist<smallest_dist)
        smallest_dist=field_dist;
        center_field_idx=i;
    end
end

% find distance from center field
for i=1:fields_num
    field_dist=euc_d(centroid_x(center_field_idx),centroid_y(center_field_idx),centroid_x(i),centroid_y(i));
    if field_dist~=0 % check that the field is not compared to itself
    if isempty(spac_exclud)==1 || isempty(find(spac_exclud==i))==1 % check index is not in exclusion list
        if use_dist_thresh==1
            if field_dist < dist_thresh % check that fields are close enough
                field_distances=[field_distances,field_dist];
            end
        else
            field_distances=[field_distances,field_dist];
        end
    end
    end
end

% find field sizes
field_sizes=[];
for i=1:fields_num
    if isempty(size_exclud)==1 || isempty(find(size_exclud==i))==1 % check index is not in exclusion list
        field_sizes=[field_sizes,length(find(fields_x(i,:)>0))];
    end
end

% find angles
angles=[];
for i=1:fields_num
    for j=1:fields_num
        angle=find_angle(centroid_x(i),centroid_y(i),centroid_x(j),centroid_y(j));
        if isnan(angle)==0
            angles=[angles,angle];
        end
    end
end
hist_binwidth=5;
angle_hist=histogram(angles,BinWidth=hist_binwidth);
%plot(angle_hist)
%disp(angle_hist.BinLimits)
%length(angle_hist.Values)
hist_vals=linspace(min(angle_hist.BinLimits),max(angle_hist.BinLimits),length(angle_hist.Values));
hist_bins=angle_hist.Values;

% find smallest angle from first centroid
cent_one_angles=[];
for i=2:length(centroid_x)
    a=find_angle(centroid_x(1),centroid_y(1),centroid_x(i),centroid_y(i));
    cent_one_angles=[cent_one_angles,a];
end

% find smallest angle from center field that is above 0 degrees
for i=1:fields_num
    if i ~= center_field_idx
    if isempty(ang_exclud)==1 || isempty(find(ang_exclud==i))==1 % check index is not in exclusion list
        a = find_angle(centroid_y(i),centroid_y(center_field_idx),centroid_x(i),centroid_x(center_field_idx));
        a = adjust_angle(a);
        if a < sml_ang_cnt_fld && a >= 0
            sml_ang_cnt_fld = a;
            sml_ang_cent_num = i;
            %fprintf("center_field_idx: %d i: %d. atan2(%f - %f, %f - %f) * 180 / pi\n",center_field_idx,i,centroid_y(i),centroid_y(center_field_idx),centroid_x(i),centroid_x(center_field_idx));
        end
    end
    end
end

mean_field_sizes=sum(field_sizes)/size(field_sizes,2);
%median_field_sizes=median(field_sizes);
mean_field_dists=sum(field_distances)/size(field_distances,2);
%median_field_dists=median(field_distances);
fprintf("Mean field sizes (area): %.2f, mean size in cm: %.2f cm\n",mean_field_sizes,mean_field_sizes*px2cm);
%fprintf("Median field sizes (area): %.2f\n",median_field_sizes);
fprintf("Mean field distances: %.2f, mean distance in cm: %.2f cm\n",mean_field_dists,mean_field_dists*px2cm);
%fprintf("Median field distances: %.2f\n",median_field_dists);
%fprintf("Minimum angle: %.2f\n",min(angles));
fprintf("Minimum angle: %.2f; Matching centroid number: %d\n",sml_ang_cnt_fld,sml_ang_cent_num);
if only_center_seven && fld_cnt~=7 fprintf("Incorrect field count for only seven fields filter.\n"); end
fprintf("Grid fields reported: %d\n",fields_num);
if only_center_seven
    fprintf("Grid fields filtered out: %d\n",filter_out);
else
    fprintf("Grid fields filtered out: %d\n",size(filter_out,2));
end
%fprintf("Grid scale score: %.2f\n",mean_field_sizes*mean_field_dists);
if print_angles
    fprintf("Frequency of angles: ");
    hist_sorted=[];
    for i=1:length(hist_vals)
        if hist_bins(i)~=0
            hist_sorted=[hist_sorted;hist_vals(i),hist_bins(i)];
        end
    end
    hist_sorted=sortrows(hist_sorted,2,'descend');
    for i=1:length(hist_sorted)
        fprintf("%.1f: %.0f; ",hist_sorted(i,1),hist_sorted(i,2));
    end
end
%{
fprintf("Centroids: ");
for i=1:length(centroid_x)
    fprintf("(%.2f,%.2f) ",centroid_x(i),centroid_y(i));
end
fprintf("\nAngles from first centroid: ");
for i=1:(length(centroid_x)-1)
    fprintf("%.2f",cent_one_angles(i));
    if i ~= (length(centroid_x)-1) fprintf(", "); end
end
%}
%{
fprintf("\nCustom angle reporting: ");
c1=1; c2=7;
a=find_angle(centroid_x(c1),centroid_y(c1),centroid_x(c2),centroid_y(c2));
fprintf("centroid %d to %d: %.2f",c1,c2,a);
%}
%fprintf("\n");

% save scores
if exist('save_field_size','var')
    if str2num([string(save_field_size)]) == 1
        fieldsize_file = fopen(filename_sizes,'at'); % append file
        % fprintf(fieldsize_file,"%f\n",mean_field_sizes);
        if only_center_seven==1
            if fld_cnt==7
	            fprintf(fieldsize_file,"%f\n",mean_field_sizes);
            else
                fprintf(fieldsize_file,"%f\n",-1); % incorrect number of fields detected
            end
        else
            fprintf(fieldsize_file,"%f\n",mean_field_sizes);
        end
	    fclose(fieldsize_file);
    end
end
if exist('save_field_spacing','var')
    if str2num([string(save_field_spacing)]) == 1
        fieldspacing_file = fopen(filename_spacings,'at'); % append file
	    % fprintf(fieldspacing_file,"%f\n",mean_field_dists);
        if only_center_seven==1
            if fld_cnt==7
	            fprintf(fieldspacing_file,"%f\n",mean_field_dists);
            else
                fprintf(fieldspacing_file,"%f\n",-1); % incorrect number of fields detected
            end
        else
            fprintf(fieldspacing_file,"%f\n",mean_field_dists);
        end
	    fclose(fieldspacing_file);
    end
end
if exist('save_field_rotation','var')
    if str2num([string(save_field_rotation)]) == 1
        fieldrotation_file = fopen(filename_rotations,'at'); % append file
	    % fprintf(fieldrotation_file,"%f\n",min(angles));
        % fprintf(fieldrotation_file,"%f\n",sml_ang_cnt_fld);
        if only_center_seven==1
            if fld_cnt==7
	            fprintf(fieldrotation_file,"%f\n",sml_ang_cnt_fld);
            else
                fprintf(fieldrotation_file,"%f\n",-1); % incorrect number of fields detected
            end
        else
            fprintf(fieldrotation_file,"%f\n",sml_ang_cnt_fld);
        end
	    fclose(fieldrotation_file);
    end
end

% plotting
size_space_ratio=(mean_field_dists/((mean_field_sizes/3.14)^0.5)/2);
if control_window_size
    if plot_fields_detected && plot_orig_firing
        f = figure;
        f.Position = [100 100 900 400];
    elseif (plot_fields_detected==1 && plot_orig_firing==0) || ...
           (plot_fields_detected==0 && plot_orig_firing==1)
        f = figure;
        f.Position = [100 100 450 400];
    end
end
if plot_fields_detected && plot_orig_firing
    tiledlayout(1,2)
    nexttile
end
if plot_orig_firing || auto_export_plots
    axis('tight')
    imagesc(heat_map_orig);
    %colormap gray;
    ylim([1 res]);
    xlim([1 res]);
    set(gca,'YDir','normal')
    if minimal_plotting_mode==0
        colorbar;
        if load_plot_from_file == 0
            plot_title=strcat('Original Autocorrelogram for Cell #',string(heat_map_selection));
        else
            plot_title=strcat('Original Autocorrelogram for Cell');
        end
        title(plot_title);
        xlabel('animal location on x axis')
        ylabel('animal location on y axis')
        set(gca,'fontsize', 14);
        axis square
    else
        if load_plot_from_file == 0
            plot_title=strcat('Cell #',string(heat_map_selection),' with Ratio=',string(size_space_ratio));
        else
            plot_title=strcat('Cell with Ratio=',string(size_space_ratio));
        end
        title(plot_title);
        set(gca,'fontsize', 14);
        axis square
    end
    if auto_export_plots==1
        ax = gca;
        if load_plot_from_file == 0
            plot_filename=strcat('../../../../holger_gridcell_theory/plots/auto_export/n',string(heat_map_selection),'_t',string(nonfld_filt_perc*100),'_rm.png');
        else
            plot_filename=strcat('../../../../holger_gridcell_theory/plots/auto_export/cell_t',string(nonfld_filt_perc*100),'_rm.png');
        end
        exportgraphics(ax,plot_filename,'Resolution',300) 
    end
end
if control_window_size && plot_fields_detected && plot_orig_firing nexttile; end
if plot_fields_detected || auto_export_plots
    gscatter(heat_map2(:,2),heat_map2(:,1),idx);
    ylim([1 res]);
    xlim([1 res]);
    set(gca,'YDir','normal')
    if minimal_plotting_mode==0
        if load_plot_from_file == 0
            plot_title=strcat('Cell #',string(heat_map_selection),' with Threshold=',string(nonfld_filt_perc),' and Ratio=',string(size_space_ratio));
        else
            plot_title=strcat('Cell with Threshold=',string(nonfld_filt_perc),' and Ratio=',string(size_space_ratio));
        end
        title(plot_title);
        xlabel('animal position on x axis')
        ylabel('animal position on y axis')
        set(gca,'fontsize', 14);
        axis square
    else
        if load_plot_from_file == 0
            plot_title=strcat('Cell #',string(heat_map_selection),' with Ratio=',string(size_space_ratio));
        else
            plot_title=strcat('Cell with Ratio=',string(size_space_ratio));
        end
        %plot_title=strcat('Cell #',string(heat_map_selection),' Field Detection Overlay');
        title(plot_title);
        set(legend,'Visible','off');
        set(gca,'fontsize', 14);
        axis square
    end
    if auto_export_plots==1
        ax = gca;
        if load_plot_from_file == 0
            plot_filename=strcat('../../../../holger_gridcell_theory/plots/auto_export/n',string(heat_map_selection),'_t',string(nonfld_filt_perc*100),'_fd.png');
        else
            plot_filename=strcat('../../../../holger_gridcell_theory/plots/auto_export/cell_t',string(nonfld_filt_perc*100),'_fd.png');
        end
        exportgraphics(ax,plot_filename,'Resolution',300) 
    end
end
if plot_legend==0 && (plot_fields_detected==1 || plot_orig_firing==1)
    if minimal_plotting_mode==0
        set(legend,'Visible','off');
    end
end

function d=euc_d(x1,y1,x2,y2)
    d=sqrt((x2-x1)^2+(y2-y1)^2); % euclidean distance
end

function a=find_angle(y2,y1,x2,x1)
    %a=atand(abs(y2-y1)/abs(x2-x1)); % angle in degrees
    a=atan2(y2 - y1, x2 - x1) * 180 / pi; % angle in degrees
end

function [filtered] = basic_convolution(image,kernel) 
    dimensions = size(image); 
    dimensions2 = size(kernel); 
 
    % define kernel center indices 
    kernelCenter_x = floor(dimensions2(1)/2); 
    kernelCenter_y = floor(dimensions2(2)/2); 
 
    image2 = zeros(dimensions(1),dimensions(2)); 
    for i = 1:dimensions(1) 
        for j = 1:dimensions(2) 
            for k = 1:dimensions2(1) 
                for l = 1:dimensions2(2) 
                    % New changes are added below 
                    ii = i+(k-kernelCenter_x); 
                    jj = j+(l-kernelCenter_y);  
                    if (ii >= 1 && ii <= dimensions(1) && jj >= 1 && jj <= dimensions(2)) 
                        %fprintf("i:%d j:%d ii:%d jj:%d k:%d l:%d\n",i,j,ii,jj,k,l);
                        image2(i,j) = image2(i,j) + image(ii,jj)* kernel(k,l); 
                    end 
                end 
            end 
        end   
    filtered = image2; 
    end 
end

function a=adjust_angle(a)
    %{
        Adjust angle to accomidate minimal angle measurement methods
        This reports a minimal angle that is not specific to the vertical
        or horizontal axis.
    %}

    if a > 45 && a <= 90
        a = 90 - a;
    elseif a > 90 && a <= 135
        a = a - 90;
    elseif a > 135 && a <= 180
        a = 180 - a;
    elseif a < -135 && a >= -180
        a = a * -1;
        a = 180 - a;
    elseif a < -90 && a >= -135
        a = a * -1;
        a = a - 90;
    elseif a < -45 && a >= -90
        a = a * -1;
        a = 90 - a;
    elseif a < 0 && a >= -45
        a = a * -1;
    end
end