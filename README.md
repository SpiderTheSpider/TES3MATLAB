# TES3MATLAB
A suite of tools for extracting and visualizing landscape data from TES3 Morrowind and its mods.  The following are one-sentence summaries of each function.

tes3matlab  
The primary code.  Use this to define options and run the other tools.

file_data = tes3matlab_extract(opts)  
Extract the data associated with the exterior landscape from the specified TES3 master or plugin file(s)


merged_data = tes3matlab_merge(file_data,opts)  
Merge the data from the multiple sources in file_data into a single data structure

merged_maps = tes3matlab_gen_maps(merged_data,opts)  
Generate full world maps from the merged data

map_img = tes3matlab_draw_maps(merged_maps,opts)  
Draw an image of a map to display as a MATLAB figure and/or export to an image file
