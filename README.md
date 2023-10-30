# Gridcell Metrics Reporting
Requirements:
<br>The DBSCAN clustering algorithm from MATLAB’s (mathworks.com) Statistics and Machine Learning Toolbox
<br>
<br>Usage:
<br>src/gridcell_metrics.m should be run in MATLAB
<br>
<br>Required Parameters:
<br>heat_map_selection: select cell to analyze (cells are numbered and file locations are stored in heat_maps_list.mat and heat_maps_ac_list.mat).
<br>
<br>Optional Parameters:
<br>nonfld_filt_perc: threshold for non-field filtering to perform. Range is 0-1.
<br>use_ac: use autocorrelogram (1) or rate map (0).
<br>only_center_seven: filter intended for autocorrelograms where fields other than the center 7 are filtered out.
<br>plot_fields_detected: plot fields the software detected.
<br>plot_orig_ratemap: plot original ratemap or autocorrelogram used as the source of data for the software.
<br>plot_legend: include legend in plot.
<br>
<br>Note: this software does not generate rate maps or autocorrelograms from animal data, although some are saved in example files with the software. Additional software is needed to create rate maps and autocorrelograms. This software uses as input the data from those plots.