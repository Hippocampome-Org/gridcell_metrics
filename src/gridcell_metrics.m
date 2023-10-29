%{
    Automatic grid cell metric reporter
    
    Author: Nate Sutton, 2023
    
    - dist_thresh: maximum field distance to be considered in group of
    field distances.
    - fld_sz_thresh: minimum field size that has contact with a border for
    the field to avoid being filtered out.
    - nonfld_filt_perc: percentage of lowest firing rates to remove from
    data. This removes non-field firing.
    
    references: https://www.mathworks.com/help/stats/dbscan.html
    https://www.mathworks.com/matlabcentral/answers/501948-how-to-automatically-find-the-angle-from-2-points
    https://www.mathworks.com/matlabcentral/answers/629408-how-to-use-convolution-on-a-2d-matrix-with-a-kernel
    http://agamtyagi.blogspot.com/2013/02/matlab-code-for-famous-mexican-hat.html
%}

load("heat_maps_list.mat"); small=35; medium=36; large=37;
load("heat_maps_ac_list.mat"); 
% real cell file_number 15 = small, 23 = medium, 17 = large

%%%%%%% parameter options %%%%%%%
heat_map_selection=large;%23;%29;%15;%23;%medium;
custom=0; custom2=0; custom3=0; custom4=0;
use_tophat_filter=0;
use_centsurr_filter=0;
% thresholds
% if heat_map_selection==small
%     dist_thresh=20; % max field distance
%     fld_sz_thresh=6; % min field neurons at border
%     nonfld_filt_perc=0.26; % low firing filter out
% elseif heat_map_selection==medium
%     dist_thresh=20;
%     fld_sz_thresh=6;
%     nonfld_filt_perc=0.26;
% elseif heat_map_selection==large
%     dist_thresh=50;
%     fld_sz_thresh=10;
%     nonfld_filt_perc=0.30;
% else
    dist_thresh=150;%15;%150;%50;%20;
    fld_sz_thresh=6;
    nonfld_filt_perc=.26;%.15;%.35;%.26;%.3;%.17;%.26;%0.26;%.3;%0.68;%0.45;
% end
dbscan_epsilon=3;
dbscan_min_pts=1;
load_custom_rm=0; % load custom rate map saved from a file
use_ac=0; % choose to use autocorrelogram (ac) (1) or standard rate map (0)
load_custom_ac=0; % load custom ac saved from a file
ac_xaxis_dim=63; % for custom ac, use this x-axis length
only_center_seven=1; % filter out fields except for seven closest to the plot center
if use_ac == 0 only_center_seven=0; end
% plotting
plot_fields_detected=1;
plot_orig_ratemap=1;
plot_legend=1; print_angles=0;
control_window_size=0;
if custom
    load_custom_rm=1;
    plot_orig_ratemap=0;
end
if custom2
    load_custom_ac=1;
    plot_orig_ratemap=0;
end
if custom3
    load_custom_rm=1;
    use_ac=0;
    only_center_seven=0;
end
if custom4
    load_custom_rm=0;
    use_ac=1;
    only_center_seven=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_fields_detected && plot_orig_ratemap control_window_size=1; end
load(heat_maps(heat_map_selection));
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
%heat_map=tophatFiltered;
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

% heat_map_orig=heat_map_flt3;
% heat_map_orig=tophatFiltered;
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
% convert to (x,y) columns
%heat_map_custom=zeros(ac_xaxis_dim);
% if use_ac
%     heat_map=fliplr(heat_map); % correct heat_map points by flipping left with right when using ac
% end
for y=1:res
    for x=1:res
        if heat_map(y,x) > 0.0
            %heat_map2=[heat_map2;[x,y]];
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
%             fields_x(i,pnt_ctr)=heat_map2(j,1);
%             fields_y(i,pnt_ctr)=heat_map2(j,2);
            fields_x(i,pnt_ctr)=heat_map2(j,2);
            fields_y(i,pnt_ctr)=heat_map2(j,1);
        end
    end
end

% filter out fields at edges without enough points
filter_out=[];
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
center_point=size(heat_map_ac,2)/2;
if only_center_seven
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
for i=1:fields_num
    for j=1:fields_num
        % check that fields are close enough
        field_dist=euc_d(centroid_x(i),centroid_y(i),centroid_x(j),centroid_y(j));
        if field_dist < dist_thresh && field_dist~=0
            field_distances=[field_distances,field_dist];
        end
    end
end

% find field sizes
field_sizes=[];
% for i=1:fields_num
%     pts_x=[];
%     pts_y=[];
%     for j=1:most_points
%         if fields_x(i,j)~=0
%             pts_x=[pts_x,fields_x(i,j)];
%             pts_y=[pts_y,fields_y(i,j)];
%         end
%     end
%     size_x=max(pts_x)-min(pts_x);
%     size_y=max(pts_y)-min(pts_y);
%     field_sizes=[field_sizes,(size_x+size_y)/2];
% end

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

mean_field_sizes=sum(field_sizes)/size(field_sizes,2);
mean_field_dists=sum(field_distances)/size(field_distances,2);
fprintf("Mean field sizes (area): %.2f\n",mean_field_sizes);
fprintf("Mean field distances: %.2f\n",mean_field_dists);
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
fprintf("Centroids ");
for i=1:length(centroid_x)
    fprintf("(%.2f,%.2f) ",centroid_x(i),centroid_y(i));
end
fprintf("\n");

% plotting
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
if plot_orig_ratemap
    axis('tight')
    imagesc(heat_map_orig);
    title('Original Rate Map');
    %colormap gray;
    ylim([1 res]);
    xlim([1 res]);
    set(gca,'YDir','normal')
    colorbar;
    xlabel('animal location on x axis')
    ylabel('animal location on y axis')
end
if control_window_size && plot_fields_detected && plot_orig_ratemap nexttile; end
if plot_fields_detected
%     gscatter(heat_map2(:,1),heat_map2(:,2),idx);
    gscatter(heat_map2(:,2),heat_map2(:,1),idx);
    title('Grid Fields Detected');
    ylim([1 res]);
    xlim([1 res]);
    set(gca,'YDir','normal')
    xlabel('neuron position on x axis')
    ylabel('neuron position on y axis')    
end
if plot_legend==0 && (plot_fields_detected==1 || plot_orig_ratemap==1)
    set(legend,'Visible','off');
end

function d=euc_d(x1,y1,x2,y2)
    d=sqrt((x2-x1)^2+(y2-y1)^2); % euclidean distance
end

function a=find_angle(x1,y1,x2,y2)
    a=atand(abs(y2-y1)/abs(x2-x1)); % angle in degrees
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