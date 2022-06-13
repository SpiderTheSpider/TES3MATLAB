function file_data = tes3matlab_extract(opts)
% TES3 MATLAB master and plugin file data extractor
% by Arin Nelson
% 06/06/2022: v0.02 init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References
% [1] https://www.mwmythicmods.com/tutorials/MorrowindESPFormat.html
%     > outdated, but is more thorough in some aspects
% [2] https://rajetic.wordpress.com/2015/03/03/terrain-methods-part-1/
%     > says clearly how to translate vertex height gradients to heights
% [3] https://en.uesp.net/wiki/Morrowind_Mod:Mod_File_Format
%     > clearer definitions but not as thorough
% [4] https://en.uesp.net/wiki/Morrowind_Mod:TES3_File_Format
%     > cites the above and possibly other useful sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIT

    % Global constants
    valid_entry_types = {'REGN','LTEX','CELL','LAND','ACTI','ALCH','APPA',...
                         'ARMO','BODY','BOOK','BSGN','CLAS','CLOT','CONT',...
                         'CREA','DIAL','DOOR','ENCH','FACT','GLOB','GMST',...
                         'INFO','INGR','LEVC','LEVI','LIGH','LOCK','MGEF',...
                         'NPC_','PGRD','PROB','RACE','REPA','SCPT','SKIL',...
                         'SNDG','SOUN','SPEL','SSCR','STAT','WEAP','MISC'};

    % Init structure containing the outputs
    file_data = struct;

    % Try to open list_file
    fid = fopen(opts.list_file,'rt');
    if(fid==-1) 
        error(['Unable to open file ' list_file ' for reading.']);
    end
    
    % Read in line by line
    file_list = textscan(fid,'%s\n','delimiter','');
    file_list = file_list{1};

    % Close file
    fclose(fid);
    clear fid;

    % Ensure files exist before attempting to read in their data
    n_file = numel(file_list);
    any_issues = false;
    for i=1:n_file

        % Attempt to open the file.
        % If successful, ensure it's a TES3 file.
        fid = fopen(file_list{i},'rb');
        if(fid<=0)
            warning(['Unable to open file ' file_list{i} ' for reading.']);
            any_issues = true;
        else
            file_signature = native2unicode( fread(fid,4,'ubit8') )';
            if(strcmp(file_signature,'TES3')~=1)
                warning(['Invalid file signature: ' file_signature ', this version only works for TES3 files!']);
                any_issues = true;
            end
            clear file_signature;
        end
        fclose(fid);

    end

    % Break code if an error was encountered and errors aren't allowed to be
    % skipped
    if(any_issues==true & opts.allow_skip==false)
        error('An issue was encountered and allow_skip was set to false.');
    end
    clear any_issues;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN

    % Loop through each file
    for i=1:n_file

        % Check if data has been read in before
        file_name = strsplit(file_list{i},{'\','/','.'});
        file_name = file_name{end-1};
        bkup_file = [opts.bkup_dir '\' file_name '.mat'];
        if(exist(bkup_file,'file')==2)
            tmp = load(bkup_file);
            file_data(i).hedr = tmp.data.hedr;
            file_data(i).regn = tmp.data.regn;
            file_data(i).ltex = tmp.data.ltex;
            file_data(i).land = tmp.data.land;
            file_data(i).cell = tmp.data.cell;
            file_data(i).file = tmp.data.file;
            clear tmp;
        else
    
            %=============================================================%
            % INIT
    
            % Open the file for reading
            fid = fopen(file_list{i},'rb');

            % Skip file signature
            nul = fread(fid,4,'ubit8');

            % Get header size
            hedr_size = fread(fid,1,'uint32');

            % Next 8 bits should be empty
            nul = fread(fid,8,'ubit8');

            % First entry should be header
            entry_name = fread(fid,4,'ubit8')';
            if( strcmp( native2unicode(entry_name), 'HEDR' )~=1 )
                error(['First entry of TES3 file should be HEDR, instead found ' native2unicode(entry_name) '!']);
            else

                % HEDR things
                hedr_size = fread(fid,1,'uint32');
                hedr_vers = fread(fid,4,'ubit8');
                hedr_type = fread(fid,4,'ubit8');   %0=esp, 1=esm, 2=ess
                hedr_comp = native2unicode(fread(fid,32,'ubit8'))';
                hedr_desc = native2unicode(fread(fid,256,'ubit8'))';
                hedr_nrec = fread(fid,1,'uint32');

                % Save header string
                file_data(i).hedr = hedr_desc;

            end
            clear hedr_size entry_name

            % Init data structures
            file_data(i).regn = struct; i_regn = 0;
            file_data(i).ltex = struct; i_ltex = 0;
            file_data(i).land = struct; i_land = 0;
            file_data(i).cell = struct; i_cell = 0;     %esp_cell(1).NAME = '';  esp_cell(1).REGION = '';  esp_cell(1).DATA = '';  esp_cell(1).COLOR = '';

            %=============================================================%
            % MAIN
            
            % While not at the end of the file, read in & parse data
            while(~feof(fid))
        
                % Get entry type and lenth
                entry_type = fread(fid,4,'ubit8')';
                if(~feof(fid))
                if(~ismember(native2unicode(entry_type),valid_entry_types))
                    warning(['Unknown entry_type: ' native2unicode(entry_type) ', scanning for valid entry_type...']);
                else

                    % Get entry info
                    entry_size  = fread(fid,1,'uint32');
                    nul         = fread(fid,4,'ubit8');    % Appears to be unused...
                    entry_flags = fread(fid,4,'ubit8');    % Unused for now

                    %-----------------------------------------------------%
                    % Switch processing method depending on entry_type
                    switch native2unicode(entry_type)

                        %.................................................% 
                        case 'REGN'                              % REGION %

                            % New region found!
                            i_regn = i_regn + 1;
                            disp(['Parsing entry REGN #' num2str(i_regn) ' (size ' num2str(entry_size) ')']);

                            % Subentries
                            mem_on = 1;
                            while mem_on < entry_size

                                % Subentry name & size
                                sub_name    = native2unicode(fread(fid,4,'ubit8')');
                                sub_size    = fread(fid,1,'uint32');
                                sub_value	= fread(fid,sub_size,'ubit8');

                                % Parse
                                switch sub_name
                                    case 'NAME';    file_data(i).regn(i_regn).name = deblank(native2unicode(sub_value)');   % ID
                                    case 'FNAM';	file_data(i).regn(i_regn).fnam = deblank(native2unicode(sub_value)');   % FULL NAME
                                    case 'CNAM';    file_data(i).regn(i_regn).cnam = sub_value(1:3);                        % COLOR (final value is alpha and is always 0, so skip)
                                    case 'SNAM';    file_data(i).regn(i_regn).snam = sub_value;                             % SOUND INFO
                                    case 'WEAT';    file_data(i).regn(i_regn).weat = sub_value;                             % WEATHER CHANCES
                                end

                                % Continue
                                mem_on = mem_on + sub_size+8; %name+size+value

                            end
                            clear mem_on sub_*;

                        %.................................................%    
                        case 'LTEX'                   % LAND TEXTURE DEFS %

                            % New land texture definition found!
                            i_ltex = i_ltex + 1;
                            disp(['Parsing entry LTEX #' num2str(i_ltex) ' (size ' num2str(entry_size) ')']);

                            % Subentries
                            mem_on = 1;
                            while mem_on < entry_size

                                % Subentry name & size
                                sub_name    = native2unicode(fread(fid,4,'ubit8')');
                                sub_size    = fread(fid,1,'uint32');

                                % Parse
                                switch sub_name
                                    case 'NAME';    file_data(i).ltex(i_ltex).name = deblank(lower(native2unicode(fread(fid,sub_size,'ubit8'))'));  % TEXTURE'S IN-GAME NAME
                                    case 'INTV';    file_data(i).ltex(i_ltex).intv = fread(fid,1,'uint32');                                         % TEXTURES'S ID IN LAND VTEX DATA
                                    case 'DATA';    file_data(i).ltex(i_ltex).data = deblank(lower(native2unicode(fread(fid,sub_size,'ubit8'))'));  % TEXTURE'S IMAGE FILE
                                end

                                % Continue
                                mem_on = mem_on + sub_size+8; %name+size+value

                            end
                            clear mem_on;

                        %.................................................%  
                        case 'CELL'                                % CELL %

                            % Memory tracker
                            mem_on = 0;

                            % Cell name
                            sub_name = native2unicode(fread(fid,4,'ubit8')');	% 'NAME'
                            sub_size = fread(fid,1,'uint32');                   % ?
                            cell_name = native2unicode(fread(fid,sub_size,'ubit8')');
                            mem_on = mem_on + 8 + sub_size;

                            % Cell data
                            sub_name  = fread(fid,4,'ubit8');    % 'DATA'
                            sub_size  = fread(fid,1,'uint32');   % 12
                            cell_data = fread(fid,4,'ubit8');
                            cell_xy   = fread(fid,2,'int32');
                            mem_on = mem_on + 20;

                            % Sometimes cells have an error
                            if(mem_on<entry_size)

                                % Region name (if it exists)
                                sub_name  = native2unicode(fread(fid,4,'ubit8')');	% 'RGNN'
                                if(strcmp(sub_name,'RGNN')==1)
                                    sub_size  = fread(fid,1,'uint32');                	% ?
                                    cell_rgnn = native2unicode(fread(fid,sub_size,'ubit8')');
                                    mem_on = mem_on + 8 + sub_size;
                                else
                                    cell_rgnn = '    ';
                                    mem_on = mem_on + 4;
                                end

                                % Determine if cell is an exterior cell by...
                                % If cell_rgnn is >4 characters
                                if(numel(cell_rgnn)>4)

                                    % Save cell info
                                    i_cell = i_cell + 1;
                                    disp(['Parsing entry CELL #' num2str(i_cell) ' (size ' num2str(entry_size) ')']);
                                    file_data(i).cell(i_cell).name = cell_name;
                                    file_data(i).cell(i_cell).rgnn = cell_rgnn;
                                    file_data(i).cell(i_cell).data = cell_data;
                                    file_data(i).cell(i_cell).xy   = cell_xy;

                                    % Continue (skip the rest for now)
                                    if(mem_on>entry_size)
                                        error('Memory record error.');   % UH OH, ERROR
                                    else
                                        nul = fread(fid,entry_size-mem_on,'ubit8');
                                    end

                                else

                                    disp(['Skipping interior CELL ' cell_name ' (size ' num2str(entry_size) ')']);
                                    if(mem_on>entry_size)
                                        error('Memory record error.');  % UH OH
                                    else
                                        nul = fread(fid,entry_size-mem_on,'ubit8');
                                    end

                                end

                            end
                            clear mem_on;

                        %.................................................%  
                        case 'LAND'                  % EXTERIOR LANDSCAPE %

                            % New land found!
                            i_land = i_land + 1;
                            disp(['Parsing entry LAND #' num2str(i_land) ' (size ' num2str(entry_size) ')']);

                            % Subentries
                            mem_on = 1;
                            while mem_on < entry_size    

                                % Subentry name & size
                                sub_name    = fread(fid,4,'ubit8')';
                                sub_size    = fread(fid,1,'uint32');

                                % Parse
                                disp([' Parsing LAND #' num2str(i_land) ' entry ' native2unicode(sub_name) '...']);
                                switch native2unicode(sub_name)

                                    % % % % % % % % % % % % % % % % % % % %
                                    case 'INTV'             % COORDINATES %          

                                        % Coordinates
                                        file_data(i).land(i_land).xy = fread(fid,2,'int32');
                                        mem_on = mem_on + 4 + 8;

                                    % % % % % % % % % % % % % % % % % % % %
                                    case 'DATA'                   % FLAGS %

                                        % Data types present
                                        nul = fread(fid,4,'ubit8');
                                        mem_on = mem_on + 4+4;

                                    % % % % % % % % % % % % % % % % % % % %
                                    case 'VNML'          % VERTEX NORMALS %
                                    if(opts.get_vnml==true)

                                        % (64x64)+1 map of RBG values, where
                                        % r = x-direction
                                        % g = y-direction
                                        % b = z-direction (vertical)
                                        % Useful for shading, computing slopes, etc.
                                        file_data(i).land(i_land).vnml = zeros(65,65,3,'int8');
                                        for ix=1:65
                                        for iy=1:65
                                        for ic=1:3
                                            file_data(i).land(i_land).vnml(ix,iy,ic) = fread(fid,1,'bit8');
                                        end
                                        end
                                        end
                                        clear ix iy ic;

                                    else

                                        % Skip
                                        file_data(i).land(i_land).vnml = zeros(65,65,3,'uint8');
                                        nul = fread(fid,12675,'bit8');

                                    end
                                    mem_on = mem_on + 4+12675;

                                    % % % % % % % % % % % % % % % % % % % %
                                    case 'VHGT'          % VECTOR NORMALS %

                                        % Offset for entire cell
                                        file_data(i).land(i_land).hoff = fread(fid,1,'float32');

                                        % 65x65 map of vertex heights
                                        file_data(i).land(i_land).vhgt = zeros(65,65,'int8');
                                        for ix=1:65
                                        for iy=1:65
                                            file_data(i).land(i_land).vhgt(ix,iy) = fread(fid,1,'int8');
                                        end
                                        end
                                        clear ix iy;

                                        % Junk
                                        nul = fread(fid,3,'ubit8');

                                        % Next
                                        mem_on = mem_on + 4+4232;

                                    % % % % % % % % % % % % % % % % % % % %
                                    case 'WNAM'                 % MINIMAP %

                                        % World minimap
                                        file_data(i).land(i_land).wnam = zeros(9,9,'uint8');
                                        for ix=1:9
                                        for iy=1:9
                                            file_data(i).land(i_land).wnam(ix,iy) = fread(fid,1,'uint8');
                                        end
                                        end
                                        clear ix iy;

                                        % Continue
                                        mem_on = mem_on + 4+81;

                                    % % % % % % % % % % % % % % % % % % % %
                                    case 'VCLR'           % VERTEX COLORS %
                                    if(opts.get_vclr==true)

                                        % 65x65 map of RBG values, where
                                        % r = x-direction
                                        % g = y-direction
                                        % b = z-direction (vertical)
                                        file_data(i).land(i_land).vclr = zeros(65,65,3,'uint8');
                                        for ix=1:65
                                        for iy=1:65
                                        for ic=1:3
                                            file_data(i).land(i_land).vclr(ix,iy,ic) = fread(fid,1,'ubit8');
                                        end
                                        end
                                        end
                                        clear ix iy ic;

                                    else

                                        % Skip
                                        file_data(i).land(i_land).vclr = zeros(65,65,3,'uint8');
                                        nul = fread(fid,12675,'bit8');

                                    end
                                    mem_on = mem_on + 4+12675;

                                    % % % % % % % % % % % % % % % % % % % %
                                    case 'VTEX'         % VERTEX TEXTURES %

                                        % 65x65 map of texture IDs
                                        file_data(i).land(i_land).vtex = zeros(16,16,'uint16');
                                        for ix=1:16
                                        for iy=1:16
                                            file_data(i).land(i_land).vtex(ix,iy) = fread(fid,1,'uint16');
                                        end
                                        end
                                        clear ix iy;

                                        % Continue
                                        mem_on = mem_on + 4+512+inf;

                                    % % % % % % % % % % % % % % % % % % % %
                                    otherwise

                                        disp([' Skipping LAND sub-entry ' sub_name ' (size ' num2str(sub_size) '...']);
                                        if(~isempty(sub_size))
                                            nul = fread(fid,sub_size,'ubit8');
                                            mem_on = mem_on + 4 + sub_size;
                                        else
                                            mem_on = mem_on + 4;
                                        end

                                end

                            end
                            clear mem_on;

                        %.................................................%
                        otherwise

                            % Skip others for now
                            if(~isempty(entry_size))
                                disp(['Skipping entry ' entry_type ' (size ' num2str(entry_size) ')...']);
                                nul = fread(fid,entry_size,'ubit8');
                            end

                    end
                    %-----------------------------------------------------%

                    % Clean-up
                    clear entry_*

                end       
                end
                
            end
            
            %=============================================================%
            % FINALE

            % When done reading file, close it
            fclose(fid);
            clear fid;

            % Oh jeez, should probably remember which file this data belongs to!
            tmp = strsplit(file_list{i},{'\','/'});
            file_data(i).file = tmp{end};
            clear tmp;

            % Backup if requested
            data = file_data(i);
            save(bkup_file,'data');
            clear data;
    
        end
    
    end
    clear i n_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end