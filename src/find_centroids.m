function [centroid_x, centroid_y]=find_centroids(com_centroids, fields_num, most_points, fields_x, fields_y, heat_map_orig)
	% find centroids
	centroid_x=[];
	centroid_y=[];
	if com_centroids == 0
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
	elseif com_centroids == 1
	    for i=1:fields_num
	        sum_x=0;
	        sum_y=0;
	        total_fr = 0;
	        for j=1:most_points
	            if fields_x(i,j)~=0
	                firing_rate = heat_map_orig(fields_y(i,j), fields_x(i,j));
	                x_pos = fields_x(i,j);
	                y_pos = fields_y(i,j);
	                total_fr = total_fr + firing_rate;
	                sum_x = sum_x + (firing_rate * x_pos);
	                sum_y = sum_y + (firing_rate * y_pos);
	            end
	        end
	        centroid_x=[centroid_x,sum_x/total_fr];
	        centroid_y=[centroid_y,sum_y/total_fr];       
	    end
	end
end