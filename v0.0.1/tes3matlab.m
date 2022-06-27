%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TES3MATLAB Master
% Extracts data from a specified list of TES3 files (ESMs and/or ESPs),
% merges the data into a single world data structure,
% generates some maps according to the user's wishes,
% then display and/or exports specified maps to image files.
% 
% by Spider the spider
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Breaker Box: Specify here which functions you want to run
% NOTE: You don't have to run these functions using this code, but it helps
%       to show the recommended order and syntax
Breaker = false(99,1);
%-------------------%
Breaker(01) = 1;    % Define options, constants, and project directories
Breaker(02) = 0;	% Load/extract data from specified list of files
Breaker(03) = 0;	% Merge data from multiple files into 1 world structure
Breaker(04) = 0;    % Generate specified world map(s)
%-------------------%
Breaker(11) = 0;    % Gather land textures
Breaker(12) = 0;    % Gather CreateMaps images
Breaker(13) = 0;    % Perform seam tests
%-------------------%
Breaker(21) = 0;    % Draw gridmap layers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS
if(Breaker(01))
    
    % Add directory with sub-functions to MATLAB'S $PATH
    %addpath('.\functions');
    
    % Init options structure
    opts = struct;

    % Project name (for data backukps and exported file naming)
    opts.project_name = 'vanilla';
    opts.project_dir  = '.\projects\vanilla';
    
    % Text file containing list of TES3 files to extract data from
    % NOTE: Make sure the load order is the same as what is used in-game!
    % NOTE: Make sure that each file name is headed by the full directory
    %       e.g. C:\games\morrowind\mods\mod_name.esp
    opts.list_file = '.\lists\filelist_vanilla.txt';
    
    % Land texture directory
    % NOTE: Currently, only file types readable by MATLAB's 'imread' are
    %       supported (bmp, jpg, png, etc.).  I first save the
    %       merged_data.ltex.data entries as a list in a text file,
    %       then use DDS Converter (available at 
    %       http://www.ddsconverter.com/) to convert these textures to PNGs
    %       after extracting them from the BSAs.
    % NOTE: In the future, this would GREATLY benefit from being able to
    %       extract data from BSAs...
    opts.ltex_dir = [opts.project_dir '\ltex'];
    
    % CreateMaps output directory
    % NOTE: This assumes that the files in list_file are what was used to
    % generate the CreateMaps output
    opts.cmap_dir = [opts.project_dir '\createmaps'];
    
    % Global options
    % NOTE: You SHOULD change these to suit your needs!  
    opts.bkup_data  = true;     % Back up extracted data so you don't have to extract it from the file every time
    opts.get_list   = {'ALL'};  % Choose 'ALL' or individual lumps (see help entry of tes3matlab_extract for more info)
    opts.map_list   = {'ALL'};  % Choose 'ALL' or individual maps (see help entry of tes3matlab_genmaps for more info)
    opts.blok_size  = 1;        % Size of a block for blockmap (in yards)
    opts.grid_size  = 40;       % Size of a cell for gridmap (in pixels)
    
    % Global constants
    % NOTE: YOU SHOULD NOT NEED TO CHANGE THESE, but if you want to adapt
    % this code for TES4 or something else, these values may change!
    opts.dims = struct;          % Contains dimension info
    opts.dims.hinc_ng = 8;       % One heightmap increment in # game units
    opts.dims.cell_ng = 8192;    % Length of cell edge in # game units
    opts.dims.cell_ny = 128;     % Length of cell edge in # yards
    opts.dims.cell_nh = 64;      % Length of a cell in heightmap grid
    opts.dims.cell_nt = 16;      % Length of a cell in land texture map grid
    opts.dims.cell_nc = 256;     % Length of a cell in a CreateMaps image file
    opts.dims.cell_nm = 9;       % Length of a cell in mini map grid
    
    % Init project directories
    opts.dir_head = opts.project_dir;           % Main project directory
    opts.dir_data = [opts.dir_head '\data'];    % Data backup directory
    opts.dir_xprt = [opts.dir_head '\exports']; % Export directory
    
    % Init project directories
    if(exist(opts.dir_head,'dir')~=7);  mkdir(opts.dir_head);   end
    if(exist(opts.dir_data,'dir')~=7);  mkdir(opts.dir_data);   end
    if(exist(opts.dir_xprt,'dir')~=7);  mkdir(opts.dir_xprt);   end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(Breaker(02));	tes3matlab_gatherdata(opts);                        end
if(Breaker(03));    merged_data = tes3matlab_mergedata(opts);           end
if(Breaker(04));	map_data    = tes3matlab_mapdata(merged_data,opts);	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if(Breaker(11));    [ltex_img, ltex_clr] = tes3matlab_get_ltex({merged_data.ltex.data},opts.ltex_dir);	end
if(Breaker(12));    [cmap_img, cmap_info] = tes3matlab_get_cmap([merged_data.land.xy],opts.cmap_dir);	end
if(Breaker(13));    [ztest_log, ctest_log, ttest_log] = tes3matlab_seamtests(merged_data);              end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tes3matlab_draw_gridmap is currently a stand-alone since it's still in development, 
% so far, maps and images have been generated individually and may be separated into separate functions in the next release
