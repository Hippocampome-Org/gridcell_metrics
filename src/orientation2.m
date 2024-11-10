function ori2_angles=orientation2(centroid_x, centroid_y, heat_map_orig, ang_exclud, fields_num)
	% find orientation angle
	angles = [];
	dists = []; % distances
	center_point = [size(heat_map_orig,1)/2,size(heat_map_orig,2)/2];
	for i=1:fields_num
	    if isempty(ang_exclud)==1 || isempty(find(ang_exclud==i))==1 % check index is not in exclusion list
	        a = find_angle(centroid_y(i),center_point(1),centroid_x(i),center_point(2)); % angle
	        a = a  * pi / 180; % convert to radians
	        %a = adjust_angle(a);
	        d = euc_d(center_point(2),center_point(1),centroid_x(i),centroid_y(i)); % distance
	        angles = [angles, a];
	        dists = [dists, d];
	    end
	end
	% find fields closest to center
	dist_indices = [1:length(dists)];
    dists_orig = dists;
	% bubble sort
	for i=1:length(dists)
		for j=1:length(dists)
			if dists(i)<dists(j)% && i ~= j
				dist_tmp = dists(i);
				dists(i)=dists(j);
				dists(j)=dist_tmp;
				dist_i_tmp = dist_indices(i);
				dist_indices(i)=dist_indices(j);
				dist_indices(j)=dist_i_tmp;
			end
		end
    end
    if length(dist_indices) <= 6
        max_cd = length(dist_indices);
    else
        max_cd = 6;
    end
    closest_dist_to_cent = dist_indices(1:max_cd);
    smallest_angles = angles(closest_dist_to_cent);
    if length(smallest_angles) <= 3
        max_sa = length(dist_indices);
    else
        max_sa = 3;
    end
    smallest_angles = smallest_angles(1:max_sa);
    smallest_angles = smallest_angles * 180 / pi; % convert to degrees
	ori2_angles = mean(smallest_angles);
end