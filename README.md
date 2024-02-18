# Gridcell Metrics Reporting

## Related article
<br>Some additional details and analyses run with this software can be found in this article.
<br>Rebecca, R. G., Ascoli, G. A., Sutton, N. M., & Dannenberg, H. (2023). Spatial periodicity in grid cell firing is explained by a neural sequence code of 2D trajectories. bioRxiv. https://doi.org/10.1101/2023.05.30.542747

## Requirement
<br>The DBSCAN clustering algorithm from MATLAB’s (mathworks.com) Statistics and Machine Learning Toolbox

## Example Usage Instructions.
<br>A user will need to supply either a rate map of grid cell firing in a 32x32 pixel matrix or an autocorrelogram in a 63x63 pixel matrix. The script ratemap2autocorr.m is supplied to convert a rate map matrix into an autocorrelogram matrix if that is wanted. The matrix should be stored in a .mat matlab file. In the matlab file, the matrix should be a matlab object named "heat_map" if it is a rate map and named "heat_map_ac" if it is an autocorrelogram.
<br>
<br>In gridcell_metrics.m, the user should enter the file path to their matrix file in the variable plot_filepath. Alternatively, a user can enter their matrix file path in the files heat_maps_list.mat or heat_maps_ac_list.mat if loading matrices from a list is preferred. The variable load_plot_from_file toggles loading a matrix from an individual file or from a list.
<br>
<br>The software includes a few required parameters and has many optional parameters. They are listed at the beginning of the script under the corresponding comment descriptions. One required parameter is px2cm, and this specifies what number of centimeters each pixel in the matrix represents, and this is used for statistics reporting. Another required parameter is nonfld_filt_perc, and this parameter sets the percentage of non-grid-field firing to filter out. Non-grid-field firing is approximated by a set percentage firing level that is below the peak firing level. For example, if the peak firing rate was 4 spikes per second, a 0.31 threshold would filter out firing that is below 4*.31=1.24 spikes per second. Suggested threshold are covered in a section below.
<br>
<br>The software can be run by clicking the run button in Matlab while the gridcell_metrics.m script is open in the editor window. Optional parameters include use_ac, and this sets if the program should process a rate map or autocorrelogram. If use_ac is set to 0 then an object in the workspace named heat_map will be used for analyses and this should be loaded from the input matrix file. If use_ac is set to 1 the program will used an object names heat_map_ac. Plotting setting include plot_fields_detected and plot_orig_firing. These should be set to show the detected field and/or the original firing plot (rate map or autocorrelogram).
<br>
<br>Statistics will be printed in the commend prompt once the program is run. These include mean field size, spacing, and angle measurements. Matching centroid number is the centroid from the detected grid field with the smallest angle to the center field. Grid fields reported is the number of grid fields included in analyses and does not take into account any manually excluded fields. Grid fields filtered out is the number of grid fields beyond the "grid fields reported" that were filtered out.

## Suggested nonfld_filt_perc thresholds
<br>...

## General Usage:
<br>src/gridcell_metrics.m should be run in MATLAB

## Required Parameters:
<br>heat_map_selection: select cell to analyze (cells are numbered and file locations are stored in heat_maps_list.mat and heat_maps_ac_list.mat).
<br>nonfld_filt_perc: threshold for non-field filtering to perform. Range is 0-1.

## Optional Parameters:
<br>use_ac: use autocorrelogram (1) or rate map (0).
<br>only_center_seven: filter intended for autocorrelograms where fields other than the center 7 are filtered out.
<br>plot_fields_detected: plot fields the software detected.
<br>plot_orig_firing: plot original ratemap or autocorrelogram used as the source of data for the software.
<br>plot_legend: include legend in plot.
<br>
<br>Note: this software does not generate rate maps from animal data, although some are saved in example files with the software. Additional software is needed to create rate maps. This software uses as input the data from those plots.