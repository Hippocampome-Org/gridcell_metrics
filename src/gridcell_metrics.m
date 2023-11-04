%{
    Automatic grid cell metric reporter
    
    Author: Nate Sutton, 2023
    
    Required parameter:
    - nonfld_filt_perc: percentage of lowest firing rates to remove from
    data. This removes non-field firing.

    Optional parameters:
    - dist_thresh: maximum field distance to be considered in group of
    field distances.
    - fld_sz_thresh: minimum field size that has contact with a border for
    the field to avoid being filtered out.
    
    references: https://www.mathworks.com/help/stats/dbscan.html
    https://www.mathworks.com/matlabcentral/answers/501948-how-to-automatically-find-the-angle-from-2-points
    https://www.mathworks.com/matlabcentral/answers/629408-how-to-use-convolution-on-a-2d-matrix-with-a-kernel
    http://agamtyagi.blogspot.com/2013/02/matlab-code-for-famous-mexican-hat.html
    https://gamedev.stackexchange.com/questions/4467/comparing-angles-and-working-out-the-difference
%}

load("heat_maps_list.mat"); small=35; medium=36; large=37;
load("heat_maps_ac_list.mat"); 
% real cell file_number 15 = small, 23 = medium, 17 = large

%%%%%%% parameter options %%%%%%%
heat_map_selection=1;%36;%23;%29;%15;%23;%medium;
if exist('heat_map_selection_custom','var')
    heat_map_selection=heat_map_selection_custom;
end
custom=0; custom2=0; custom3=0; custom4=0;
use_tophat_filter=0;
use_centsurr_filter=0;
minimal_plotting_mode=0;
if exist('minimal_plotting_mode_custom','var')
    minimal_plotting_mode=minimal_plotting_mode_custom;
end
auto_export_plots=0;
if exist('auto_export_plots_custom','var')
    auto_export_plots=auto_export_plots_custom;
end
filename_sizes='saved_results/field_size_records.txt';
if exist('filename_sizes_custom','var')
    filename_sizes=filename_sizes_custom;
end
filename_spacings='saved_results/field_spacing_records.txt';
if exist('filename_spacings_custom','var')
    filename_spacings=filename_spacings_custom;
end
filename_rotations='saved_results/save_field_rotation.txt';
if exist('filename_rotations_custom','var')
    filename_rotations=filename_rotations_custom;
end
% threshold
nonfld_filt_perc=.26;
if exist('saved_nonfld_filt_perc','var')
    nonfld_filt_perc=saved_nonfld_filt_perc;
end
% optional settings
use_dist_thresh=0;
dist_thresh=500;
use_fld_sz_thresh=0;
fld_sz_thresh=6;
% end
dbscan_epsilon=3;
dbscan_min_pts=1;
load_custom_rm=0; % load custom rate map saved from a file
use_ac=1; % choose to use autocorrelogram (ac) (1) or standard rate map (0)
load_custom_ac=0; % load custom ac saved from a file
ac_xaxis_dim=63; % for custom ac, use this x-axis length
only_center_seven=1; % filter out fields except for seven closest to the plot center
if use_ac == 0 only_center_seven=0; end
% plotting
plot_fields_detected=1;
plot_orig_ratemap=1;
if exist('plot_fields_detected_custom','var') plot_fields_detected=plot_fields_detected_custom; end
if exist('plot_orig_ratemap_custom','var') plot_orig_ratemap=plot_orig_ratemap_custom; end
plot_legend=0; print_angles=0;
control_window_size=0;
if custom load_custom_rm=1; plot_orig_ratemap=0; end
if custom2 load_custom_ac=1; plot_orig_ratemap=0; end
if custom3 load_custom_rm=1; use_ac=0; only_center_seven=0; end
if custom4 load_custom_rm=0; use_ac=1; only_center_seven=1; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_fields_detected && plot_orig_ratemap control_window_size=1; end
if exist('use_newly_generated_rm','var')
    % only load heat map from file if no heat map is already in memory
    if str2num([string(use_newly_generated_rm)]) == 0
        load(heat_maps(heat_map_selection));
    end
else
    load(heat_maps(heat_map_selection));
end
if load_custom_rm
    load("heat_maps_real_custom_list.mat");
    load(heat_maps_real_custom_list(heat_map_selection));
    heat_map=heat_map_custom;
end
if use_ac
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
        non_zero=find(fields_x(i,:)~=0);
        non_zero_sz=size(non_zero,2);
        fld_x=fields_x(non_zero);
        fld_y=fields_y(non_zero);
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
    closest_seven=closest_seven(1:7);
    % filter out unwanted fields
    filter_out=(fields_num-7);
    fields_num=7;
    fields_x=fields_x(closest_seven,:);
    fields_y=fields_y(closest_seven,:);
    centroid_x=centroid_x(:,closest_seven);
    centroid_y=centroid_y(:,closest_seven);
    % only keep points in the seven fields
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
        if use_dist_thresh==1
            if field_dist < dist_thresh % check that fields are close enough
                field_distances=[field_distances,field_dist];
            end
        else
            field_distances=[field_distances,field_dist];
        end
    end
end

%{
% commented out because development was not completed on this
% find surrounding field spacing
% for i=1:fields_num
%     if i ~= center_field_idx % avoid comparisons with center field because those have already been done
%         smallest_dist=euc_d(centroid_x(i),centroid_y(i),centroid_x(1),centroid_y(1));
%         smallest_dist_idx=1;
%         for j=1:fields_num
%             field_dist=euc_d(centroid_x(i),centroid_y(i),centroid_x(j),centroid_y(j));
%             if field_dist~=0 && field_dist<smallest_dist
%                 smallest_dist=field_dist;
%                 smallest_dist_idx=j;
%             end
%         end
%         field_distances=[field_distances,smallest_dist];
%     end
% end

% find fields with similar angles to center field
fields_similar_angles=[];
fields_similar_angles_values=[];
fields_similar_angles_values_2=[];
for i=1:fields_num
    if i ~= center_field_idx
    angle_field_1=find_angle(centroid_x(i),centroid_y(i),centroid_x(center_field_idx),centroid_y(center_field_idx));
    field_2_idx=i+1; if field_2_idx>fields_num field_2_idx=field_2_idx-fields_num; end
    if field_2_idx==center_field_idx field_2_idx=field_2_idx+1; end
    if field_2_idx>fields_num field_2_idx=field_2_idx-fields_num; end
    angle_field_2=find_angle(centroid_x(field_2_idx),centroid_y(field_2_idx),centroid_x(center_field_idx),centroid_y(center_field_idx));
    % normalize angle difference
    norm_angle_difference = 180 - abs(abs(angle_field_1 - angle_field_2) - 180);
    smallest_dist=norm_angle_difference;
    smallest_dist_idx=field_2_idx;
    for j=1:fields_num
        if j ~= center_field_idx
        angle_field_1=find_angle(centroid_x(i),centroid_y(i),centroid_x(center_field_idx),centroid_y(center_field_idx));
        angle_field_2=find_angle(centroid_x(j),centroid_y(j),centroid_x(center_field_idx),centroid_y(center_field_idx));
        norm_angle_difference = 180 - abs(abs(angle_field_1 - angle_field_2) - 180);
        if (i ~= j && norm_angle_difference < smallest_dist)
            % check to avoid adding fields pair already included
            skip=0;
            if size(fields_similar_angles,1)>0
                for k=1:size(fields_similar_angles,1)
                    f1=fields_similar_angles(k,1);
                    f2=fields_similar_angles(k,2);
                    if (i==f1 && j==f2) || (i==f2 && j==f1)
                        skip=1;
                    end
                end
            end
            if skip==0
                smallest_dist = norm_angle_difference;
                smallest_dist_idx = j;
            end
        end
%       commented out because development is incomplete
%         % direction of angles
%         norm_angle_1_direction=angle_field_1;
%         norm_angle_2_direction=angle_field_2;
%         % this assumes a difference no more than 90 degrees with angle 2 when 
%         % angle 1 is between 0 and 90 degrees
%         if (angle_field_1 < 90 && angle_field_2 > 270) norm_angle_2_direction=norm_angle_2_direction-360;
%         if (angle_field_2 < 90 && angle_field_1 > 270) norm_angle_1_direction=norm_angle_2_direction-360;
%         direction_sign=norm_angle_2_direction-norm_angle_1_direction;
%         direction=1; % counter-clockwise direction
%         if direction_sign>0
%             direction=0; % clockwise direction
%         end
        end
    end
    fields_similar_angles=[fields_similar_angles; [i,smallest_dist_idx]];
    % report angles
    angle_field_1=find_angle(centroid_x(i),centroid_y(i),centroid_x(center_field_idx),centroid_y(center_field_idx));
    angle_field_2=find_angle(centroid_x(smallest_dist_idx),centroid_y(smallest_dist_idx),centroid_x(center_field_idx),centroid_y(center_field_idx));
    norm_angle_difference = 180 - abs(abs(angle_field_1 - angle_field_2) - 180);
    fields_similar_angles_values=[fields_similar_angles_values, norm_angle_difference];
    fields_similar_angles_values_2=[fields_similar_angles_values_2;[angle_field_1,angle_field_2]];
    % % %
    field_dist=euc_d(centroid_x(i),centroid_y(i),centroid_x(smallest_dist_idx),centroid_y(smallest_dist_idx));
    field_distances=[field_distances,field_dist];
    end
end
%}

% find field sizes
field_sizes=[];
for i=1:fields_num
    field_sizes=[field_sizes,length(find(fields_x(i,:)>0))];
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

mean_field_sizes=sum(field_sizes)/size(field_sizes,2);
%median_field_sizes=median(field_sizes);
mean_field_dists=sum(field_distances)/size(field_distances,2);
%median_field_dists=median(field_distances);
fprintf("Mean field sizes (area): %.2f\n",mean_field_sizes);
%fprintf("Median field sizes (area): %.2f\n",median_field_sizes);
fprintf("Mean field distances: %.2f\n",mean_field_dists);
%fprintf("Median field distances: %.2f\n",median_field_dists);
fprintf("Minimum angle: %.2f\n",min(angles));
fprintf("Grid fields reported: %d\n",fields_num);
if only_center_seven
    fprintf("Grid fields filtered out: %d\n",filter_out);
else
    fprintf("Grid fields filtered out: %d\n",size(filter_out,2));
end
fprintf("Grid scale score: %.2f\n",mean_field_sizes*mean_field_dists);
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
fprintf("Centroids: ");
for i=1:length(centroid_x)
    fprintf("(%.2f,%.2f) ",centroid_x(i),centroid_y(i));
end
fprintf("\nAngles from first centroid: ");
for i=1:(length(centroid_x)-1)
    fprintf("%.2f",cent_one_angles(i));
    if i ~= (length(centroid_x)-1) fprintf(", "); end
end
%{
fprintf("\nCustom angle reporting: ");
c1=1; c2=7;
a=find_angle(centroid_x(c1),centroid_y(c1),centroid_x(c2),centroid_y(c2));
fprintf("centroid %d to %d: %.2f",c1,c2,a);
%}
fprintf("\n");

% save scores
if exist('save_field_size','var')
    if str2num([string(save_field_size)]) == 1
        fieldsize_file = fopen(filename_sizes,'at'); % append file
	    fprintf(fieldsize_file,"%f\n",mean_field_sizes);
	    fclose(fieldsize_file);
    end
end
if exist('save_field_spacing','var')
    if str2num([string(save_field_spacing)]) == 1
        fieldspacing_file = fopen(filename_spacings,'at'); % append file
	    fprintf(fieldspacing_file,"%f\n",mean_field_dists);
	    fclose(fieldspacing_file);
    end
end
if exist('save_field_rotation','var')
    if str2num([string(save_field_rotation)]) == 1
        fieldrotation_file = fopen(filename_rotations,'at'); % append file
	    fprintf(fieldrotation_file,"%f\n",min(angles));
	    fclose(fieldrotation_file);
    end
end

% plotting
size_space_ratio=(mean_field_dists/((mean_field_sizes/3.14)^0.5)/2);
if control_window_size
    if plot_fields_detected && plot_orig_ratemap
        f = figure;
        f.Position = [100 100 900 400];
    elseif (plot_fields_detected==1 && plot_orig_ratemap==0) || ...
           (plot_fields_detected==0 && plot_orig_ratemap==1)
        f = figure;
        f.Position = [100 100 450 400];
    end
end
if plot_fields_detected && plot_orig_ratemap
    tiledlayout(1,2)
    nexttile
end
if plot_orig_ratemap || auto_export_plots
    axis('tight')
    imagesc(heat_map_orig);
    %colormap gray;
    ylim([1 res]);
    xlim([1 res]);
    set(gca,'YDir','normal')
    if minimal_plotting_mode==0
        colorbar;
        plot_title=strcat('Original Rate Map for Cell #',string(heat_map_selection));
        title(plot_title);
        xlabel('animal location on x axis')
        ylabel('animal location on y axis')
        set(gca,'fontsize', 14);
        axis square
    else
        plot_title=strcat('Cell #',string(heat_map_selection),' with Ratio=',string(size_space_ratio));
        title(plot_title);
        set(gca,'fontsize', 14);
        axis square
    end
    if auto_export_plots==1
        ax = gca;
        plot_filename=strcat('../../../../holger_gridcell_theory/plots/auto_export/n',string(heat_map_selection),'_t',string(nonfld_filt_perc*100),'_rm.png');
        exportgraphics(ax,plot_filename,'Resolution',300) 
    end
end
if control_window_size && plot_fields_detected && plot_orig_ratemap nexttile; end
if plot_fields_detected || auto_export_plots
    gscatter(heat_map2(:,2),heat_map2(:,1),idx);
    ylim([1 res]);
    xlim([1 res]);
    set(gca,'YDir','normal')
    if minimal_plotting_mode==0
        plot_title=strcat('Cell #',string(heat_map_selection),' with Threshold=',string(nonfld_filt_perc),' and Ratio=',string(size_space_ratio));
        title(plot_title);
        xlabel('animal position on x axis')
        ylabel('animal position on y axis')
        set(gca,'fontsize', 14);
        axis square
    else
        plot_title=strcat('Cell #',string(heat_map_selection),' with Ratio=',string(size_space_ratio));
        %plot_title=strcat('Cell #',string(heat_map_selection),' Field Detection Overlay');
        title(plot_title);
        set(legend,'Visible','off');
        set(gca,'fontsize', 14);
        axis square
    end
    if auto_export_plots==1
        ax = gca;
        plot_filename=strcat('../../../../holger_gridcell_theory/plots/auto_export/n',string(heat_map_selection),'_t',string(nonfld_filt_perc*100),'_fd.png');
        exportgraphics(ax,plot_filename,'Resolution',300) 
    end
end
if plot_legend==0 && (plot_fields_detected==1 || plot_orig_ratemap==1)
    if minimal_plotting_mode==0
        set(legend,'Visible','off');
    end
end

function d=euc_d(x1,y1,x2,y2)
    d=sqrt((x2-x1)^2+(y2-y1)^2); % euclidean distance
end

function a=find_angle(x1,y1,x2,y2)
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