function [gs_orientation, gs_orientations_std]=extract_grid_orientation(orientations)
	orientations = mod(orientations,60);
    corr_orientations = [];
    for i=1:length(orientations)
        % For every angle extract min to 60 deg
        diff_60 = orientations(i) - 60;
        if abs(diff_60) < abs(orientations(i))
            corr_orientations = [corr_orientations, diff_60];
        else
            corr_orientations = [corr_orientations, orientations(i)];
        end
    end

    % Check 30 degree flips 
    if mean(abs(abs(corr_orientations) - 30)) < mean(abs(corr_orientations))
        % Yes, angles close to 30 degrees (flipping axis)
        % Try to reach consensus
        if median(corr_orientations) < 0 % Make everything negative
        	for i=1:length(corr_orientations)
        		if corr_orientations(i) > 0 corr_orientations(i) = corr_orientations(i) * -1; end
        	end
        else % Make everything positive
        	for i=1:length(corr_orientations)
        		if corr_orientations(i) < 0 corr_orientations(i) = corr_orientations(i) * -1; end
        	end
        end
    end

    % Extract average and standard deviation
    orientation     = mean(corr_orientations);
    orientation_std = std(corr_orientations);

    % Test / correct orientation 
    test_matrix = [abs(orientation-60); abs(orientation)];
    % find the index of the minimum value in the test_matrix. Check if that index is 1.
    if find(test_matrix==min(test_matrix)) == 1
    	orientation = orientation - 60;
    end

    gs_orientation = orientation;
    gs_orientations_std = orientation_std;
end