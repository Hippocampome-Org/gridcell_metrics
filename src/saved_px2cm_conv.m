% saved pixel to centimeter conversion settings

if load_plot_from_file==0
if heat_map_selection < 30
    px2cm=100/32;
elseif heat_map_selection >= 30 && heat_map_selection <= 34
    px2cm=45/32;
elseif heat_map_selection >= 35 && heat_map_selection <= 100
    px2cm=100/32;
end
end