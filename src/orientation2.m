function ori2_angles=orientation2(centroid_x, centroid_y, heat_map_orig, ang_exclud, fields_num, min_orientation2)
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

	% where two fields have a very similar orientation, discard the more distant one
	orient_distsq = abs(circ_dist2(angles));
	close_fields = orient_distsq(find(orient_distsq < min_orientation2));
	close_fields = triu(close_fields, 1); % Upper triangle only - set lower triangle to zero
                                          % k=1: +1 offset from diagonal: set diagonal to zero too
    to_del = [];
    for i=1:size(close_fields,1)
    	for i=2:size(close_fields,2)
    		if dists(i) > dists(j)
    			to_del = [to_del, i];
    		else
    			to_del = [to_del, j];
    		end
    	end
    end
    orig_indices = linspace(1,length(dists),length(dists));
    to_del_inv = setdiff(orig_indices,to_del);
    dists = dists(to_del_inv);
    angles = angles(to_del_inv);

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

    % save at most the 6 closest fields
    if length(dist_indices) <= 6
        max_cd = length(dist_indices);
    else
        max_cd = 6;
    end
    closest_dist_to_cent = dist_indices(1:max_cd);
    
    % save at most the 3 smallest angles among the 6 closest fields
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