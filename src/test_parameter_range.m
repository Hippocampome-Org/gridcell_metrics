save_field_size=1;
save_field_spacing=1;

for i=20:33
    fprintf("\nnow processing %d%%\n",i);
    filename_sizes_custom=strcat('saved_results/field_size_records_',string(i),'perc.txt');
    filename_spacings_custom=strcat('saved_results/field_spacings_records_',string(i),'perc.txt');
    saved_nonfld_filt_perc=i*.01;
    for j=1:29
        fprintf("\nnow processing cell %d\n",j);
        heat_map_selection_custom=j;
        gridcell_metrics;
    end
end