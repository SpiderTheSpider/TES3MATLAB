%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TES3MATLAB Master
% Extracts data from a specified list of TES3 files (ESMs and/or ESPs),
% generate some exterior landscape maps according to the user's wishes,
% then display and/or export maps to image files.
% 
% by Arin Nelson
% 06/04/2022: v0.01 init
% 06/06/2022: v0.02 init
% 06/13/2022: v0.03 init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Breaker Box: Specify here which functions you want to run
% NOTE: You don't have to run these functions using this code, but it helps
%       to show the recommended call syntax and order
Breaker = false(99,1);
%-------------------%
Breaker(01) = 0;	% Load/extract data from specified list of files
Breaker(02) = 1;	% Merge data from multiple files into 1 data structure
Breaker(03) = 0;    % Generate specified world map(s)
Breaker(04) = 0;    % Draw specified map(s)
Breaker(05) = 0;    % Display and/or export specified map(s)
%-------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
    
    % First things first, add dir with sub-functions to MATLAB'S $PATH
    addpath('.\functions');
    
    % Options structure
    opts = struct;
    
    % Project name (for backup and exported file naming)
    opts.project_name = 'vanilla';
    
    % Text file containing list of TES3 files to extract data from
    % NOTE: Make sure the load order is the same as what is used in-game!
    % NOTE: Make sure that each file name is headed by the full directory
    %       e.g. C:\games\morrowind\mods\mod_name.esp
    opts.list_file = '.\lists\filelist_vanilla.txt';
    
    % Directory to look for land texture files
    % NOTE: This, and all other defined directories, need to be manually
    %       created!  This code will not create them for security reasons.
    % NOTE: Currently, only file types readable by MATLAB's 'imread' are
    %       supported (bmp, jpg, png, etc.).  I use tes3matlab_export_ltex
    %       to generate a list of land textures used by the imported
    %       file(s), then use DDS Converter (available at 
    %       http://www.ddsconverter.com/) to convert these textures to PNGs
    opts.ltex_dir = '.\images\ltex_pngs\';
    
    % CreateMaps output directory
    % NOTE: This assumes that the same files used in the CreateMaps output
    %       are what is listed in input_file
    opts.createmaps_dir = '.\images\createmaps\vanilla\';
    
    % GridMap file
    % NOTE: Unused at the moment, but could be another source of info...
    % NOTE: Open and save the master file as a .TIFF so MATLAB can read it!
    opts.gridmap_file = '.\images\gridmap\Tamriel_MasterImage_v061122.tiff';

    % Export directory, where results are saved to
    opts.export_dir = '.\exports\vanilla\';
    
    % Backup directory, where extracted data is stored so you don't have to
    % run the extraction every time
    opts.bkup_dir  = '.\backups\';
    
    % Global options
    % NOTE: You SHOULD change these to suit your needs!  
    % NOTE: It GREATLY speeds things up if you pick-and-choose which datums
    %       to extract, process, and display!
    opts.allow_skip = false;    % Allow skipping the extraction for a data file if an error is encountered
    opts.bkup_data  = true;     % Back up extracted data so you don't have to extract it from the file every time
    opts.get_vclr	= false;	% Extract vertex color map of LAND entries
    opts.get_vnml	= true; 	% Extract vertex normal map of LAND entries
    opts.map_to_gen = {'ALL'};  % Choose 'ALL' or individual maps (see help entry of tes3matlab_gen_maps for more info)
    opts.dims.cell_bgrid = 16;  % Used to generate a 'block map', a map where dx=dy=dz (think Minecraft)
    
    % Global constants
    % NOTE: YOU SHOULD NOT NEED TO CHANGE THESE, but if you want to adapt
    % this code for TES4 or something else, these values may change!
    opts.dims = struct;             % Contains the various dimension definitions
    opts.dims.cell_gu    = 8192;	% Length of cell edge in game units
    opts.dims.cell_yd    = 128;     % Length of cell edge in yards
    opts.dims.h_inc      = 8;       % One heightmap increment (in vertical) is this many game units
    opts.dims.cell_hgrid = 64;      % Length of a cell in heightmap grid
    opts.dims.cell_tgrid = 16;      % Length of a cell in land texture map grid
    opts.dims.cell_cgrid = 16;      % Length of a cell in CreateMaps image file
    opts.dims.cell_ggrid = 40;      % Length of a cell in GridMap image file
    opts.dims.cell_mgrid = 9;       % Length of a cell in mini map grid

% END OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(01)); file_data = tes3matlab_extract(opts);                  end
if(Breaker(02)); merged_data = tes3matlab_merge(file_data,opts);        end
if(Breaker(03)); map_data = tes3matlab_genmaps(merged_data,opts);       end
if(Breaker(04)); map_imgs = tes3matlab_drawmaps(map_data,opts);         end
if(Breaker(05)); tes3matlab_exportmaps(map_imgs,opts);                  end
