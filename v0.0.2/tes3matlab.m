%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TES3MATLAB Master
% Extracts data from a specified list of TES3 files (ESMs and/or ESPs),
% merges the data into a single world data structure,
% generates world maps according to the specified options,
% then exports maps as layers in a .tif file.
% 
% Note:
% Individual files are expected to be in 32-bit architecture (meaning max
% file size is 4Gb, variables expected to be at max uint32 or float/single,
% etc.), but combined/merged data should be able to go beyond that.
% 
% by Spider the spider
% 06/04/2022: v0.0.0 init
% 06/06/2022: v0.0.1 init
% 06/13/2022: v0.0.2 init
% 06/24/2022: v0.0.3 init
% 06/27/2022: v0.0.4 init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Breaker Box: Specify here which functions you want to run
% NOTE: You don't have to run these functions using this code, but it helps
%       to show the recommended order and syntax
Breaker = false(99,1);
%-------------------%
Breaker(01) = 1;    % Definitions and initializations
Breaker(02) = 1;    % Gather data
%-------------------%
Breaker(03) = 0;    % Clean gathered data
Breaker(04) = 0;    % Backup gathered data
%-------------------%
Breaker(05) = 1;	% Merge data from multiple files into 1 data structure
Breaker(06) = 1;    % Map data onto world grid (includes seam testing)
%-------------------%
Breaker(07) = 1;    % Draw raw maps        
Breaker(08) = 0;    % Draw automap layers   
%-------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(01))

    %=====================================================================%
    % Init structures
    % NOTE: Don't need to init variables in MATLAB like you do in cpp,
    % etc., but I do anyways to make converting to another code easier
    
    % Init options structure
    opts = struct;
    opts.proj_name   = '';              % User-defined project name
    opts.proj_dir    = '';              % Directory to store project files
    opts.proj_infile = '';              % Text list of input files
    opts.ltex_list   = '';              % Text list of land texture colors
    opts.ltex_dir    = '';              % Directory of land texture images
    opts.ltex_size   = uint16(0);       % Size of landscape texture images
    opts.cmap_dir    = '';              % Directory of CreateMaps images
    opts.cmap_size   = uint16(0);       % Size of CreateMaps images
    opts.gmap_size   = uint16(0);       % Gridmap scale in pixels
    opts.gmap_lclr   = uint8([0 0 0]);  % Gridmap land color
    opts.gmap_wclr   = uint8([0 0 0]);  % Gridmap water color
    opts.bkup_raw    = false;           % Export (backup) extracted data
    opts.bkup_merged = false;           % Export (backup) merged data
    opts.bkup_mapped = false;           % Export (backup) mapped data
    opts.bkup_imaged = false;           % Export (backup) image data
    opts.do_cleaning = false;           % Clean extracted data
    opts.data_dir    = '';
    opts.xprt_dir    = '';

    % Init constants structure
    cnst         = struct;      % Contains dimensions & conversion factors
    cnst.hinc_ng = uint16(0);	% Heightmap increment in game units
    cnst.cell_ng = uint16(0);	% Length of cell edge in game units
    cnst.cell_ny = uint16(0);	% Length of cell edge in yards
    cnst.cell_nh = uint16(0);	% Length of cell in heightmap grid
    cnst.cell_nt = uint16(0);	% Length of cell in landtex grid
    cnst.cell_nc = uint16(0);	% Length of cell in CreateMaps image
    cnst.cell_nm = uint16(0);	% Length of cell in minimap grid
     
    %=====================================================================%
    % Define options
    
    % Project name (for data backukps and exported file naming)
    opts.proj_name = 'mapdump';
    opts.proj_dir  = '.\projects\mapdump';

    % Text file containing list of TES3 files to extract data from
    % NOTE: Make sure the load order is the same as what is used in-game!
    % NOTE: Make sure that each file name is headed by the full directory
    %       e.g. C:\games\morrowind\mods\mod_name.esp
    %       (case does not matter, all is converted to lowercase anyways)
    opts.proj_infile = '.\lists\filelist_mapdump.txt';
    
    % List of colors to use for textures in ltex
    % NOTE: In the interest of time, I've precompiled a list of  colors 
    % for land textures used by vanilla, Tamriel Rebuilt, and Project
    % Tamriel, computed from the mean color of each image.  Since MATLAB
    % cannot natively read BSA archives or DDS images, I would need a lot
    % of time to do this operation programatically, hence why I developed
    % my own list independently.  For more details, see the ReadMe.
    opts.lclr_file = '.\lists\ltex_colors_raw.txt';
    opts.ltex_dir  = 'F:\OSOM_Data_Repo\git\_transfer\TR\Tools\tes3matlab\images\ltex_bmp';
    opts.ltex_size = uint16(16);    % Square size of a ltex image
    
    % CreateMaps output directory
    % NOTE: To get the CreateMaps output, you must load the files as listed
    % in list_file in vanilla Morrowind, then run this command in the
    % in-game console as soon as possible (before Juib even speaks):
    % CreateMaps "Morrowind.esm"
    opts.cmap_dir     = 'F:\OSOM_Data_Repo\git\_transfer\TR\Tools\tes3matlab\v0.0.3\projects\tamriel_rebuilt\createmaps';
    opts.cmap_size(1) = 256;   % Square size of a CreateMaps image
    
    % Backup options
    opts.bkup_raw    = true;	% Backup data after extracting
    opts.bkup_merged = true;   % Backup data after merging
    opts.bkup_mapped = true;   % Backup data after mapping
    opts.bkup_imaged = false;  % Backup data after imaging
    
    % Gridmap options
    opts.gmap_size(1) = 40;     % Length of cell side in Gridmap in pixels
    opts.gmap_lclr(:) = [222 199 176];  % Gridmap land color (RGB)
    opts.gmap_wclr(:) = [133 164 179];  % Gridmap water color (RGB)
    
    % Other options
    %opts.file_to_show = [4];       % Only draw maps for these files
    opts.draw_xr      = [-29 -9];   % x cell range for drawing maps
    opts.draw_yr      = [ 14 34];   % y cell range for drawing maps
    opts.do_cleaning  = true;       % Clean data after it is extracted?
    %opts.merge_land   = logical([0 0 0 0 1]); % Include land data in merge?
    opts.merge_land   = true(43,1);
    opts.rclr_file    = '.\lists\road_colors.txt';  
    
    % Options for displaying cell seams
    opts.seam_grid  = 8;    % One cell side will be this # pixels
    opts.seam_width = 1;	% Thickness of line at cell border with seam
    %opts.seam_color = [63 0 0; 255 192 192];    % Colorscale for seam size
    
    % Derived values
    opts.data_dir = [opts.proj_dir '\data'];    % Data backup directory
    opts.xprt_dir = [opts.proj_dir '\exports']; % Image export directory
    
    %=====================================================================%
    % Define constants   
    
    % Global constants
    % NOTE: YOU SHOULD NOT NEED TO CHANGE THESE, but if you want to adapt
    % this code for TES4 or something else, these values may change!
    cnst.hinc_ng = 8;       % Heightmap increment in game units
    cnst.cell_ng = 8192;	% Length of cell edge in game units
    cnst.cell_ny = 128;     % Length of cell edge in# yards
    cnst.cell_nh = 64;      % Length of cell in heightmap grid
    cnst.cell_nt = 16;      % Length of cell in land tex grid
    cnst.cell_nc = 256;     % Length of cell in CreateMaps img
    cnst.cell_nm = 9;       % Length of cell in mini map grid
    cnst.n_weat  = 10;      % Possible # of weather entries

    %=====================================================================%
    % Other stuff
    
    % Add path containing secondary functions
    addpath('.\functions');
    
    % Add paths containing tertiary functions
    addpath('.\functions\tes3matlab_parserec');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(02)==true)
    
    % Initialize project directories
    if(exist(opts.proj_dir,'dir')~=7);  mkdir(opts.proj_dir);   end
    if(exist(opts.data_dir,'dir')~=7);  mkdir(opts.data_dir);   end
    if(exist(opts.xprt_dir,'dir')~=7);  mkdir(opts.xprt_dir);   end
    
    % Gather raw data
    file_data = tes3matlab_gatherdata(opts,cnst);  
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(03)==true)
    file_data = tes3matlab_clean(file_data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(04)==true)
    
    % Loop through files
    for i=1:numel(file_data)
        
        % Get this file's name and use it to generate the data file name
        this_file = strsplit(file_data(i).file,{'\','/','.'});
        this_file = this_file{end-1};
        save_file = [opts.data_dir '\' this_file '.mat'];   % NOTE: IN FUTURE WILL BE NETCDF FOR MORE UNIVERSAL ACCESS!
        
        % Save data
        tes3matlab_savedata(file_data(i),save_file); %,opts); 
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(05)==true)
    save_file = [opts.data_dir '\_merged.mat'];
    if(exist(save_file,'file')==2)
        merged_data = tes3matlab_loaddata(save_file,{'regn','cell','ltex','land'});
    else
        merged_data = tes3matlab_mergedata(file_data,opts,cnst);
        tes3matlab_savedata(merged_data,save_file);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(06)==true)
    save_file = [opts.data_dir '\_mapped.mat'];
    if(exist(save_file,'file')==2)
        mapped_data = tes3matlab_loaddata(save_file,{'xr','yr','avail'...
                                          'regn','name','vtex','vhgt',...
                                          'vnml','vclr','fhgt','fclr',...
                                          'wnam','seams'});
    else
        mapped_data = tes3matlab_mapdata(file_data,merged_data,opts,cnst);
        tes3matlab_savedata(mapped_data,save_file);
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(07)==true)
    %ltex = tes3matlab_get_ltex({merged_data.ltex.data},opts.ltex_dir,opts.lclr_file);
    imaged_data = tes3matlab_drawdata(mapped_data,merged_data,opts,cnst);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(08)==true)
    %cmap = tes3matlab_get_cmap([merged_data.land.x;merged_data.land.y],opts.cmap_dir);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
