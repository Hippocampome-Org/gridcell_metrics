%{
    This file includes which grid fields were manually selected for 
    exclusion for the purposes of statistical analyses of cells.
%}

% based on nonfld_filt_perc=.31;
if load_plot_from_file==0
if heat_map_selection==3
    spac_exclud = [6,7]; 
    size_exclud = [6,7];
    ang_exclud = [6,7];
elseif heat_map_selection==17
    spac_exclud = []; 
    size_exclud = [5,7];
    ang_exclud = [];
elseif heat_map_selection==29
    spac_exclud = [6,7]; 
    size_exclud = [6,7];
    ang_exclud = [];
elseif heat_map_selection==30
    spac_exclud = [3,4,5,6,7];%[3,4,6,7]; 
    size_exclud = [2,3,4,5,6,7];%[3,4,6,7];
    ang_exclud = [];
elseif heat_map_selection==31
    spac_exclud = [2,3,6,7]; 
    size_exclud = [2,3,4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==32
    spac_exclud = [];%[4,5,6,7]; 
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==37
    spac_exclud = []; 
    size_exclud = [2,3,4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==39
    spac_exclud = [4,6];
    size_exclud = [4,5,6,7];
    ang_exclud = [4,6];
elseif heat_map_selection==43
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==44
    spac_exclud = [4,6];
    size_exclud = [4,5,6,7];
    ang_exclud = [4,6];
elseif heat_map_selection==45
    spac_exclud = [];
    size_exclud = [2,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==46
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==47
    spac_exclud = [];
    size_exclud = [2,3,4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==48
    spac_exclud = [];
    size_exclud = [2,3,4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==49
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==50
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==51
    spac_exclud = [];
    size_exclud = [2,3,4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==52
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==53
    spac_exclud = [];
    size_exclud = [6,7];
    ang_exclud = [];
elseif heat_map_selection==54
    spac_exclud = [];
    size_exclud = [6,7];
    ang_exclud = [];
elseif heat_map_selection==78
    spac_exclud = [];
    size_exclud = [2,3,4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==79
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==80
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==81
    spac_exclud = [];
    size_exclud = [2,3,4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==86
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
elseif heat_map_selection==89
    spac_exclud = [];
    size_exclud = [4,5,6,7];
    ang_exclud = [];
end
end