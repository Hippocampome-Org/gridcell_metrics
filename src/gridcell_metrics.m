%{
    Automatic grid cell metric reporter
    
    Author: Nate Sutton, 2023
    
    references: https://www.mathworks.com/help/stats/dbscan.html
    https://www.mathworks.com/matlabcentral/answers/501948-how-to-automatically-find-the-angle-from-2-points
    https://www.mathworks.com/matlabcentral/answers/629408-how-to-use-convolution-on-a-2d-matrix-with-a-kernel
    http://agamtyagi.blogspot.com/2013/02/matlab-code-for-famous-mexican-hat.html
    https://gamedev.stackexchange.com/questions/4467/comparing-angles-and-working-out-the-difference
    https://stackoverflow.com/questions/12801400/find-the-center-of-mass-of-points
%}

function gridcell_metrics(output_filename,px2cm,nonfld_filt_perc,load_plot_from_file,plot_filepath,use_binary_input,heat_map_selection,use_ac,convert_to_ac,use_dist_thresh,dist_thresh,use_fld_sz_thresh,fld_sz_thresh,manual_field_exclud,load_px2cm_conv,dbscan_epsilon,dbscan_min_pts,load_custom_rm,load_custom_ac,only_center_seven,only_center_seven_inout_excl,use_tophat_filter,use_centsurr_filter,minimal_plotting_mode,auto_export_plots,filename_sizes,filename_spacings,filename_rotations,plot_fields_detected,plot_orig_firing,plot_legend,print_angles,control_window_size,custom,custom2,custom3,custom4,sml_ang_cnt_fld,sml_ang_cent_num,advanced_detection,advanced_detection_maxdist,advanced_detection_ang_inc,com_centroids,in_out_fields_ratio,cov_between_fields_reporting,report_gridscore,report_orientation2,min_orientation2,spac_exclud,size_exclud,ang_exclud)

px2cm = str2num(px2cm);
nonfld_filt_perc = str2num(nonfld_filt_perc);
load_plot_from_file = str2num(load_plot_from_file);
use_binary_input = str2num(use_binary_input);
heat_map_selection = str2num(heat_map_selection);
use_ac = str2num(use_ac);
convert_to_ac = str2num(convert_to_ac);
use_dist_thresh = str2num(use_dist_thresh);
dist_thresh = str2num(dist_thresh);
use_fld_sz_thresh = str2num(use_fld_sz_thresh);
fld_sz_thresh = str2num(fld_sz_thresh);
manual_field_exclud = str2num(manual_field_exclud);
load_px2cm_conv = str2num(load_px2cm_conv);
dbscan_epsilon = str2num(dbscan_epsilon);
dbscan_min_pts = str2num(dbscan_min_pts);
load_custom_rm = str2num(load_custom_rm);
load_custom_ac = str2num(load_custom_ac);
only_center_seven = str2num(only_center_seven);
only_center_seven_inout_excl = str2num(only_center_seven_inout_excl);
use_tophat_filter = str2num(use_tophat_filter);
use_centsurr_filter = str2num(use_centsurr_filter);
minimal_plotting_mode = str2num(minimal_plotting_mode);
auto_export_plots = str2num(auto_export_plots);
plot_fields_detected = str2num(plot_fields_detected);
plot_orig_firing = str2num(plot_orig_firing);
plot_legend = str2num(plot_legend);
print_angles = str2num(print_angles);
control_window_size = str2num(control_window_size);
custom = str2num(custom);
custom2 = str2num(custom2);
custom3 = str2num(custom3);
custom4 = str2num(custom4);
sml_ang_cnt_fld = str2num(sml_ang_cnt_fld);
sml_ang_cent_num = str2num(sml_ang_cent_num);
advanced_detection = str2num(advanced_detection);
advanced_detection_maxdist = str2num(advanced_detection_maxdist); % max distance to test
advanced_detection_ang_inc = str2num(advanced_detection_ang_inc); % angle to increment
com_centroids = str2num(com_centroids); % use center of mass to find centroids
in_out_fields_ratio = str2num(in_out_fields_ratio); % find ratio of firing intensity within fields compared to outside of them
cov_between_fields_reporting = str2num(cov_between_fields_reporting); % report standard deviation of mean firing rates between fields
report_gridscore = str2num(report_gridscore);
report_orientation2 = str2num(report_orientation2);
min_orientation2 = str2num(min_orientation2);
spac_exclud = str2double(strsplit(spac_exclud,"|"));
size_exclud = str2double(strsplit(size_exclud,"|"));
ang_exclud = str2double(strsplit(ang_exclud,"|"));

%disp(heat_map_selection)
disp("Please wait while results are being generated. This may take several minutes.")

%% Load configuration settings %%
%config;
% config_filename = "config.txt";
% spac_exclud = []; size_exclud = []; ang_exclud = []; % specify any fields to exclude from statistics
% config_file = fopen(config_filename);
% config_vars = textscan(config_file,'%s',3,'Delimiter','|');
% config_data = textscan(config_file,'%s %s %s',37,'Delimiter',',');
% fclose(config_file);
% for i=1:length(config_data{1})
%     var_name = string(config_data{1}(i));        
%     var_data = string(config_data{2}(i));
%     command = sprintf("%s = %s;",var_name,var_data);
%     eval(command);
% end
%%%
load("heat_maps_list.mat"); small=35; medium=36; large=37; % list of rate map plots. real cell file_number 15 = small, 23 = medium, 17 = large
load("heat_maps_ac_list.mat"); % list of autocorrelogram plots
if load_px2cm_conv==1 saved_px2cm_conv; end
if exist('heat_map_selection_custom','var') heat_map_selection=heat_map_selection_custom; end
if use_ac == 0 only_center_seven=0; end
if exist('minimal_plotting_mode_custom','var') minimal_plotting_mode=minimal_plotting_mode_custom; end
if exist('auto_export_plots_custom','var') auto_export_plots=auto_export_plots_custom; end
if exist('filename_sizes_custom','var') filename_sizes=filename_sizes_custom; end
if exist('filename_spacings_custom','var') filename_spacings=filename_spacings_custom; end
if exist('filename_rotations_custom','var') filename_rotations=filename_rotations_custom; end
if exist('saved_nonfld_filt_perc','var') nonfld_filt_perc=saved_nonfld_filt_perc; end
if exist('plot_fields_detected_custom','var') plot_fields_detected=plot_fields_detected_custom; end
if exist('plot_orig_firing_custom','var') plot_orig_firing=plot_orig_firing_custom; end
if custom load_custom_rm=1; plot_orig_firing=0; end
if custom2 load_custom_ac=1; plot_orig_firing=0; end
if custom3 load_custom_rm=1; use_ac=0; only_center_seven=0; end
if custom4 load_custom_rm=0; use_ac=1; only_center_seven=1; end
%%%

if manual_field_exclud==1
    excluded_fields_filename = sprintf("excluded_fields/heatmap_%d.txt", heat_map_selection);
    excluded_fields_file = fopen(excluded_fields_filename);
    excluded_fields_vars = textscan(excluded_fields_file,'%s',3,'Delimiter','|');
    excluded_fields_data = textscan(excluded_fields_file,'%d %d %d');
    fclose(excluded_fields_file);
    for i=1:length(excluded_fields_vars{1})
        var_data = "[";
        for j=1:length(excluded_fields_data{i})
            var_data = sprintf("%s%s",var_data,string(excluded_fields_data{i}(j)));
            if j == (length(excluded_fields_data{i})-1)
                var_data = sprintf("%s,",var_data);
            else
                var_data = sprintf("%s]",var_data);
            end
        end
        var_name=string(excluded_fields_vars{1}(i));
        command = sprintf("%s = %s;",var_name,var_data);
        eval(command);
    end
    %excluded_fields; 
end

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
    if use_binary_input == 0
        heat_map=readmatrix(plot_filepath);
    else if use_binary_input == 1
        load(plot_filepath);
    end
    if exist('heat_map_ac','var') heat_map=heat_map_ac; end
end
if use_ac 
    if exist('heat_map_ac','var') 
        heat_map=heat_map_ac; 
    end 
end
if convert_to_ac
    use_ac = 1;
    heat_map_ac = moserac(heat_map);
    heat_map=heat_map_ac;
end
if use_ac && load_plot_from_file==0
    load(heat_maps_ac_list(heat_map_selection));
    heat_map=heat_map_ac;
end
if load_custom_ac
    load("heat_maps_ac_custom_list.mat");
    load(heat_maps_ac_custom_list(str2num(heat_map_selection)));
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
    heat_map_nonfilt = heat_map;
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
[centroid_x, centroid_y]=find_centroids(com_centroids, fields_num, most_points, fields_x, fields_y, heat_map_orig);

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
    if only_center_seven_inout_excl == 1;
        all_fields=linspace(1,fields_num,fields_num);
        not_closest_seven=setdiff(all_fields,closest_seven);
        not_closest_seven_x=[];
        not_closest_seven_y=[];
        for i=1:fields_num
            non_closest_found = 0;
            for j=1:length(not_closest_seven)
                if i == not_closest_seven(j)
                    non_closest_found = 1;
                end
            end
            if non_closest_found == 1
                not_closest_seven_x = [not_closest_seven_x, fields_x(i,:)];
                not_closest_seven_y = [not_closest_seven_y, fields_y(i,:)];
            end
        end
    end
        %fields_x_orig=fields_x;
        %fields_y_orig=fields_y;

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

% store field points for plotting
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

% advanced and localized detection of fields
if advanced_detection == 1
    test=[];
    centroid_fr = []; % centroid firing rates
    fields_x = [];
    fields_y = [];
    % find firing rates of centroids
    for i=1:length(centroid_x)
        centroid_fr = [centroid_fr,heat_map_nonfilt(round(centroid_y(i)),round(centroid_x(i)))];
        %fprintf("x: %d; y: %d; r: %f\n",round(centroid_x(i)),round(centroid_y(i)),centroid_fr(i));
    end
    % detect boundary area of fields. scan 360 degrees at 5 degree increments. check firing rate
    % starting at center in degree direction using 1 pixel x and y increments. furthest distance
    % is set to advanced_detection_maxdist
    x_inside = []; % x point inside boundary
    y_inside = []; % y point inside boundary
    y_max = size(heat_map_nonfilt,1);
    x_max = size(heat_map_nonfilt,2);
    for ci=1:size(centroid_x,2)
        xc = centroid_x(ci); % center x
        yc = centroid_y(ci); % center y
        x_inside_temp = []; % x point inside boundary
        y_inside_temp = []; % y point inside boundary
        for a=0:advanced_detection_ang_inc:360
            x_prior = xc;
            y_prior = yc;
            boundary_found = false;
            for h=0:advanced_detection_maxdist
                [x_scan,y_scan]=find_ver_hor(a, h);
                x_scan = round(x_scan + xc);
                y_scan = round(y_scan + yc);
                if x_scan > x_max x_scan = y_max; end
                if y_scan > y_max y_scan = y_max; end
                if x_scan < 1 x_scan = 1; end
                if y_scan < 1 y_scan = 1; end
                fr_scan = heat_map_nonfilt(y_scan,x_scan);
                if fr_scan < (centroid_fr * nonfld_filt_perc)
                    if boundary_found == false
                        x_inside_temp = [x_inside_temp,x_prior];
                        y_inside_temp = [y_inside_temp,y_prior];
                        boundary_found = true;
                    end
                end
                x_prior = x_scan;
                y_prior = y_scan;
            end
        end
        % fill in field area based on borders
        % scan on y-axis
        y_boundary_max = max(y_inside_temp);
        y_boundary_min = min(y_inside_temp);
        x_inside_temp2 = [];
        y_inside_temp2 = [];
        for by=y_boundary_min:y_boundary_max
            test_indices = find(y_inside_temp==by);
            min_v = min(x_inside_temp(test_indices));
            max_v = max(x_inside_temp(test_indices));
            for bx=min_v:max_v
                x_inside_temp2 = [x_inside_temp2, bx];
                y_inside_temp2 = [y_inside_temp2, by];
                x_inside = [x_inside, bx];
                y_inside = [y_inside, by];
            end
        end
        % scan on x-axis
        x_boundary_max = max(x_inside_temp);
        x_boundary_min = min(x_inside_temp);
        for bx=x_boundary_min:x_boundary_max
            test_indices = find(x_inside_temp==bx);
            min_v = min(y_inside_temp(test_indices));
            max_v = max(y_inside_temp(test_indices));
            for by=min_v:max_v
                x_inside_temp2 = [x_inside_temp2, bx];
                y_inside_temp2 = [y_inside_temp2, by];
                x_inside = [x_inside, bx];
                y_inside = [y_inside, by];
            end
        end
        % remove duplicates
        x_inside_temp3 = [];
        y_inside_temp3 = [];
        for di=1:length(x_inside_temp2)
            x1 = x_inside_temp2(di); y1 = y_inside_temp2(di);
            match_counter = 0;
            duplicate_found = false;
            for di2=1:length(x_inside_temp3)
                x2 = x_inside_temp3(di2); y2 = y_inside_temp3(di2);
                if x1 == x2 && y1 == y2
                    match_counter = match_counter + 1;
                end
                if match_counter > 0
                    duplicate_found = true;
                end
            end
            if duplicate_found == false
                x_inside_temp3 = [x_inside_temp3, x_inside_temp2(di)];
                y_inside_temp3 = [y_inside_temp3, y_inside_temp2(di)];
            end
        end
        % x_inside_temp3 = x_inside_temp;
        % y_inside_temp3 = y_inside_temp;
        for ii=1:length(x_inside_temp3)
            fields_x(ci,ii) = x_inside_temp3(ii);
            fields_y(ci,ii) = y_inside_temp3(ii);
            if ci == 1
                test = [test; x_inside_temp3(ii),y_inside_temp3(ii)];
            end
        end
    end

    % store field points for plotting
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

    % filter out non-field area
    heat_map=zeros(size(heat_map_nonfilt,1));
    for i=1:size(x_inside,2)
            xi = x_inside(i);
            yi = y_inside(i);
            heat_map(yi,xi)=heat_map_nonfilt(yi,xi);
        end
    end

    % create filtered heat map for saving to a file
    heat_map_custom=heat_map;
    heat_map_custom(find(heat_map_custom>0))=1.1111;
    heat_map_custom(find(heat_map_custom<0))=0;
    heat_map=heat_map_custom;

    % update values
    fields_num=max(idx);
    points_per_field = histcounts(idx);
    most_points=max(points_per_field);

    % update centroids
    [centroid_x, centroid_y]=find_centroids(com_centroids, fields_num, most_points, fields_x, fields_y, heat_map_orig);
end

% find field spacing
field_distances=[];
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

% find firing intensity ratio within fields compared to outside of them
if in_out_fields_ratio == 1
    out_field_firing = [];
    in_field_firing = [];
    for x=1:size(heat_map_orig,2)
        for y=1:size(heat_map_orig,1)
            % standard within field search
            in_field_found = false;
            for i=1:size(fields_x,1)
                for j=1:length(fields_x(i,:))
                    if x == fields_x(i,j) && y == fields_y(i,j) && isnan(heat_map_orig(y,x)) == 0
                        in_field_firing = [in_field_firing, heat_map_orig(y,x)];
                        in_field_found = true;
                    end
                end
            end
            % field search excluding fields outside of central seven if that option is enabled
            not_closest_seven_found = false;
            if only_center_seven == 1 && only_center_seven_inout_excl == 1
                for i=1:size(not_closest_seven_x,2)
                    if x == not_closest_seven_x(i) && y == not_closest_seven_y(i) ...
                        && isnan(heat_map_orig(y,x)) == 0
                        not_closest_seven_found = true;
                    end
                end
            end
            if in_field_found == false && not_closest_seven_found == false && isnan(heat_map_orig(y,x)) == 0
                out_field_firing = [out_field_firing, heat_map_orig(y,x)];
            end
        end
    end
    mean_in_field_firing = mean(in_field_firing);
    mean_out_field_firing = mean(out_field_firing);
    in_out_fields_ratio_value = mean_in_field_firing / mean_out_field_firing;
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
if print_angles
    hist_binwidth=5;
    angle_hist=histogram(angles,BinWidth=hist_binwidth);
    %plot(angle_hist)
    %disp(angle_hist.BinLimits)
    %length(angle_hist.Values)
    hist_vals=linspace(min(angle_hist.BinLimits),max(angle_hist.BinLimits),length(angle_hist.Values));
    hist_bins=angle_hist.Values;
    title("Histogram of Angles");
    xlabel('Angle Between Fields')
    ylabel('Number of Angles')
    set(gca,'fontsize', 14);
    axis square
end

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

% alternative orientation reporting
if report_orientation2 == 1
    ori2_angles=orientation2(centroid_x, centroid_y, heat_map_orig, ang_exclud, fields_num, min_orientation2);
    [gs_orientation, gs_orientations_std]=extract_grid_orientation(ori2_angles);
end

% evaluate the firing rate per field
if cov_between_fields_reporting == 1
    cov_frs = [];
    for i=1:size(fields_x,1)
        non_zero_indices = find(fields_x(i,:)~=0);
        frs = [];
        for j=1:length(non_zero_indices)
            fr = heat_map_orig(fields_y(i,j),fields_x(i,j));
            frs = [frs, fr];
        end
        mean_fr = mean(frs);
        std_fr = std(frs);
        cov_frs = [cov_frs, mean_fr/std_fr];
    end
    cov_mean = mean(cov_frs);
end

% grid score
if report_gridscore == 1
    m = [];
    [HDgridScore,gridness3Score]=get_HDGridScore(m,m,m,heat_map);
end

% report metrics
output_file = fopen(output_filename,'w');
mean_field_sizes=sum(field_sizes)/size(field_sizes,2);
%median_field_sizes=median(field_sizes);
mean_field_dists=sum(field_distances)/size(field_distances,2);
%median_field_dists=median(field_distances);
output_string = sprintf("metric,value\nmean_field_sizes,%.2f\nmean_field_sizes_px2cm,%.2f\n",mean_field_sizes,mean_field_sizes*px2cm);
fprintf(output_file,output_string);
output_string = sprintf("Mean field sizes (area): %.2f; mean size in cm: %.2f cm\n",mean_field_sizes,mean_field_sizes*px2cm);
fprintf(output_string);
%fprintf("Median field sizes (area): %.2f\n",median_field_sizes);
output_string = sprintf("mean_field_dists,%.2f\nmean_field_dists_px2cm,%.2f\n",mean_field_dists,mean_field_dists*px2cm);
fprintf(output_file,output_string);
output_string = sprintf("Mean field distances: %.2f; mean distance in cm: %.2f cm\n",mean_field_dists,mean_field_dists*px2cm);
fprintf(output_string);
%fprintf("Median field distances: %.2f\n",median_field_dists);
%fprintf("Minimum angle: %.2f\n",min(angles));
output_string = sprintf("minimum angle,%.2f\nmatching_centroid_number,%d\n",sml_ang_cnt_fld,sml_ang_cent_num);
fprintf(output_file,output_string);
output_string = sprintf("Minimum angle: %.2f; Matching centroid number: %d\n",sml_ang_cnt_fld,sml_ang_cent_num);
fprintf(output_string);
if only_center_seven && fld_cnt~=7 
    output_string = sprintf("Error: incorrect field count for only seven fields filter.\n");
    fprintf(output_file,output_string); fprintf(output_string); 
end
output_string = sprintf("grid_fields_reported,%d\n",fields_num);
fprintf(output_file,output_string); 
output_string = sprintf("Grid fields reported: %d; ",fields_num);
fprintf(output_string);
if only_center_seven
    output_string = sprintf("grid_fields_filtered_out,%d\n",filter_out);
    fprintf(output_file,output_string); 
    output_string = sprintf("Grid fields filtered out: %d\n",filter_out);
    fprintf(output_string);
else
    output_string = sprintf("grid_fields_filtered_out,%d\n",size(filter_out,2));
    fprintf(output_file,output_string); 
    output_string = sprintf("Grid fields filtered out: %d\n",size(filter_out,2));
    fprintf(output_string);
end
%size_space_ratio=(mean_field_dists/((mean_field_sizes/3.14)^0.5)/2);
size_space_ratio=(mean_field_dists/(((mean_field_sizes/3.14)^0.5)*2));
output_string = sprintf("size_space_ratio,%.4f\n",size_space_ratio);
fprintf(output_file,output_string);
output_string = sprintf("Spacing and Size Ratio: %.2f\n",size_space_ratio);
fprintf(output_string);
if in_out_fields_ratio == 1
    output_string = sprintf("Inside and outside field firing intensity ratio: %.2f\n",in_out_fields_ratio_value);
    fprintf(output_string);
    output_string = sprintf("in_out_ratio,%.4f\n",in_out_fields_ratio_value);
    fprintf(output_file,output_string);
end
if cov_between_fields_reporting == 1
    output_string = sprintf("Mean Coefficient of Variation Between Grid Fields: %.2f\n",cov_mean);
    fprintf(output_string);
    output_string = sprintf("cov_mean,%.4f\n",cov_mean);
    fprintf(output_file,output_string);    
end
if report_gridscore == 1
    output_string = sprintf("HD grid score: %.2f; Gridness3Score: %.2f\n",HDgridScore,gridness3Score);
    fprintf(output_string);
    output_string = sprintf("hd_grid_score,%.4f\ngridness_3_score,%.4f\n",HDgridScore,gridness3Score);
    fprintf(output_file,output_string); 
end
if report_orientation2 == 1
    output_string = sprintf("Orientation2 angle: %.2f; standard deviation of angles: %.2f\n",gs_orientation,gs_orientations_std);
    fprintf(output_string);
    output_string = sprintf("gs_orientation,%.4f\ngs_orientations_std,%.4f\n",gs_orientation,gs_orientations_std);
    fprintf(output_file,output_string); 
end
%fprintf("Grid scale score: %.2f\n",mean_field_sizes*mean_field_dists);
if print_angles 
    output_string = sprintf("Frequency of angles: ");
    fprintf(output_string);
    hist_sorted=[];
    for i=1:length(hist_vals)
        if hist_bins(i)~=0
            hist_sorted=[hist_sorted;hist_vals(i),hist_bins(i)];
        end
    end
    hist_sorted=sortrows(hist_sorted,2,'descend');
    for i=1:length(hist_sorted)
        output_string = sprintf("frequency_of_angle,%.1f,%.0f\n; ",hist_sorted(i,1),hist_sorted(i,2));
        fprintf(output_file,output_string); 
        output_string = sprintf("%.1f: %.0f; ",hist_sorted(i,1),hist_sorted(i,2));
        fprintf(output_string);
    end
    fprintf("\n");
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
fclose(output_file);

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
            plot_filename=strcat('auto_export/n',string(heat_map_selection),'_t',string(nonfld_filt_perc*100),'_rm.png');
        else
            plot_filename=strcat('auto_export/cell_t',string(nonfld_filt_perc*100),'_rm.png');
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
            plot_filename=strcat('auto_export/n',string(heat_map_selection),'_t',string(nonfld_filt_perc*100),'_fd.png');
        else
            plot_filename=strcat('auto_export/cell_t',string(nonfld_filt_perc*100),'_fd.png');
        end
        exportgraphics(ax,plot_filename,'Resolution',300) 
    end
end
if plot_legend==0 && (plot_fields_detected==1 || plot_orig_firing==1)
    if minimal_plotting_mode==0
        set(legend,'Visible','off');
    end
end

exitcode = 0;

end