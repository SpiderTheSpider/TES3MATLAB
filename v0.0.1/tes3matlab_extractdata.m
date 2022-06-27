function tes3matlab_extractdata(input_file,output_file)
% tes3matlab_extractdata(input_file,output_file,opts)
% Extracts exterior landmass data from TES3 master (.esm) or plugin (.esp)
% file in_file and saves it to the NetCDF file out_file.
% 
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
                         'SNDG','SOUN','SPEL','SSCR','STAT','WEAP','MISC',...
                         'MAST','DATA'};

    % Init data structures
    data      = struct;
    data.hedr = struct;
	data.regn = struct;  i_regn = 0;
	data.ltex = struct;  i_ltex = 0;
	data.land = struct;  i_land = 0;
	data.cell = struct;  i_cell = 0;     

	% Open the file for reading
	fid = fopen(input_file,'rb');

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
        data.hedr.size = fread(fid,1,'uint32');
        data.hedr.vers = fread(fid,4,'ubit8');
        data.hedr.type = fread(fid,4,'ubit8');   %0=esp, 1=esm, 2=ess
        data.hedr.comp = native2unicode(fread(fid,32,'ubit8'))';
        data.hedr.desc = native2unicode(fread(fid,256,'ubit8'))';
        data.hedr.nrec = fread(fid,1,'uint32');

	end
    clear hedr_size entry_name

    % Init log
    disp(['Extracting data from file ' input_file '...']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN LOOP
    entry_type = fread(fid,4,'ubit8')';
    while(~feof(fid))
        
        % Continue if entry_type is valid
        % If it's a MAST entry, its format is a little different from the
        % others, so parse it separately.
        if(~ismember(native2unicode(entry_type),valid_entry_types))
            warning([' Unknown entry_type: ' native2unicode(entry_type) ', scanning for valid entry_type...']);
        elseif(strcmp(native2unicode(entry_type),'MAST')==1)
                        
            % Get master name
            entry_size  = fread(fid,1,'uint32'); 
            nul         = fread(fid,entry_size,'ubit8');
                        
            % Get master info
            nul = fread(fid,4,'ubit8');     % String that says DATA
            nul = fread(fid,4,'ubit8');     % Size of next entry (always 8)
            nul = fread(fid,8,'ubit8');     % Size of master file
        
        else    
            
            % Information on this entry
            entry_size  = fread(fid,1,'uint32');    
            nul         = fread(fid,4,'ubit8');    % Appears to be unused...
            entry_flags = fread(fid,4,'ubit8');    % Unused for now
        
            %=============================================================%
            % Switch processing method depending on entry_type
            switch native2unicode(entry_type)
                
                %.........................................................% 
                case 'REGN'                                      % REGION %

                    % New region found!
                    i_regn = i_regn + 1;
                    disp([' Parsing entry REGN #' num2str(i_regn) ' (size ' num2str(entry_size) ')']);

                    % Subentries
                    mem_on = 1;
                    while mem_on < entry_size

                        % Subentry name & size
                        sub_name  = native2unicode(fread(fid,4,'ubit8')');
                        sub_size  = fread(fid,1,'uint32');
                        mem_on    = mem_on + 8;
                        
                        % Parse
                        switch sub_name
                        	case 'NAME';	data.regn(i_regn).name = deblank(native2unicode(fread(fid,sub_size,'ubit8'))');   % ID
                            case 'FNAM';	data.regn(i_regn).fnam = deblank(native2unicode(fread(fid,sub_size,'ubit8'))');   % FULL NAME
                            case 'CNAM';    data.regn(i_regn).cnam = fread(fid,sub_size,'ubit8');	% COLOR (final value is alpha and is always 0, so skip)
                            case 'SNAM';    data.regn(i_regn).snam = fread(fid,sub_size,'ubit8');	% SOUND INFO
                            case 'WEAT';    data.regn(i_regn).weat = fread(fid,sub_size,'ubit8');	% WEATHER CHANCES
                            case 'BNAM';    nul = fread(fid,sub_size,'ubit8');                      % ID OF CREATURE THAT YOU'LL ENCOUNTER IF INTERRUPTED WHILE RESTING
                            otherwise;
                                error(['  Unknown REGN sub-entry ' sub_name]);
                        end

                        % Continue
                        mem_on = mem_on + sub_size; %name+size+value

                    end
                    clear mem_on sub_*;
                
                %.........................................................%    
                case 'LTEX'                               % LAND TEXTURES %

                	% New land texture definition found!
                    i_ltex = i_ltex + 1;
                    disp([' Parsing entry LTEX #' num2str(i_ltex) ' (size ' num2str(entry_size) ')']);

                    % Switch though possible entries and skip DELE
                    do_dele = false;
                    mem_on  = 1;
                    while mem_on < entry_size

                    	% Subentry name & size
                        sub_name = native2unicode(fread(fid,4,'ubit8')');
                        sub_size = fread(fid,1,'uint32');

                        % Parse
                        switch sub_name
                        	case 'NAME';    data.ltex(i_ltex).name = deblank(lower(native2unicode(fread(fid,sub_size,'ubit8'))'));  % TEXTURE'S IN-GAME NAME
                            case 'INTV';    data.ltex(i_ltex).intv = fread(fid,1,'uint32');                                         % TEXTURES'S ID IN LAND VTEX DATA
                            case 'DATA';    data.ltex(i_ltex).data = deblank(lower(native2unicode(fread(fid,sub_size,'ubit8'))'));  % TEXTURE'S IMAGE FILE
                            case 'DELE';    do_dele = true;     nul = fread(fid,sub_size,'ubit8');
                            otherwise;      
                                error(['Unexpected LTEX sub-entry: ' sub_name]);
                        end

                        % Continue
                        mem_on = mem_on + 8+sub_size; %name+size+value

                    end
                    clear mem_on sub_*;

                    % Delete this entry if specified
                    if(do_dele)
                    	i_ltex = i_ltex - 1;
                        data.ltex = data.ltex(1:end-1);
                    end
                    clear do_dele;
                     
                %.........................................................%  
                case 'CELL'                                        % CELL %

                    % Memory tracker
                    mem_on = 0;

                    % Cell name
                    sub_name  = native2unicode(fread(fid,4,'ubit8')');	% 'NAME'
                    sub_size  = fread(fid,1,'uint32');                   % ?
                    cell_name = native2unicode(fread(fid,sub_size,'ubit8')');
                    mem_on    = mem_on + 8 + sub_size;

                    % Cell data
                    sub_name  = fread(fid,4,'ubit8');    % 'DATA'
                    sub_size  = fread(fid,1,'uint32');   % 12
                    cell_data = fread(fid,4,'ubit8');
                    cell_xy   = fread(fid,2,'int32');
                    mem_on    = mem_on + 20;

                    % Sometimes cells have an error
                    if(mem_on<entry_size)

                        % Region name (if it exists)
                        sub_name  = native2unicode(fread(fid,4,'ubit8')');	% 'RGNN'
                        if(strcmp(sub_name,'RGNN')==1)
                            sub_size  = fread(fid,1,'uint32');                	% ?
                            cell_rgnn = native2unicode(fread(fid,sub_size,'ubit8')');
                            mem_on    = mem_on + 8 + sub_size;
                        else
                            cell_rgnn = '    ';
                            mem_on    = mem_on + 4;
                        end

                        % Determine if cell is an exterior cell by...
                        % If cell_rgnn is >4 characters
                        if(numel(cell_rgnn)>4)

                        	% Save cell info
                            i_cell = i_cell + 1;
                            disp([' Parsing entry CELL #' num2str(i_cell) ' (size ' num2str(entry_size) ')']);
                            data.cell(i_cell).name = cell_name;
                            data.cell(i_cell).rgnn = cell_rgnn;
                            data.cell(i_cell).data = cell_data;
                            data.cell(i_cell).xy   = cell_xy;

                            % Continue (skip the rest for now)
                            nul = fread(fid,entry_size-mem_on,'ubit8');

                        else

                            % Skip data for interior cell
                         	disp(['  Skipping, is interior cell...']);
                            nul = fread(fid,entry_size-mem_on,'ubit8');

                        end

                    end
                    clear mem_on sub_* cell_*;               
                    
                %.........................................................%  
                case 'LAND'                          % EXTERIOR LANDSCAPE %

                	% New land found!
                    i_land = i_land + 1;
                    disp([' Parsing entry LAND #' num2str(i_land) ' (size ' num2str(entry_size) ')']);

                    % Subentries
                    mem_on = 1;
                    while mem_on < entry_size    

                    	% Subentry name & size
                        sub_name    = fread(fid,4,'ubit8')';
                        sub_size    = fread(fid,1,'uint32');
                        mem_on      = mem_on + 8;

                        % Parse
                        disp(['  Parsing LAND #' num2str(i_land) ' entry ' native2unicode(sub_name) ' (size ' num2str(sub_size) ')...']);
                        switch native2unicode(sub_name)

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            case 'INTV'                     % COORDINATES %          

                            	% Coordinates
                                data.land(i_land).xy = fread(fid,2,'int32');
                                mem_on = mem_on + 8;

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            case 'DATA'                           % FLAGS %

                            	% Data types present
                                nul = fread(fid,1,'uint32');
                                mem_on = mem_on + 4;

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            case 'VNML'                  % VERTEX NORMALS %
                                
                            	% (64x64)+1 map of RBG values, where
                                % r = x-direction
                                % g = y-direction
                                % b = z-direction (vertical)
                                % Useful for shading, computing slopes, etc.
                                data.land(i_land).vnml = zeros(65,65,3,'int8');
                                for ix=1:65
                                for iy=1:65
                                	data.land(i_land).vnml(ix,iy,:) = fread(fid,3,'bit8');
                                end
                                end
                                clear ix iy;

                                % Continue
                                mem_on = mem_on + 12675;

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            case 'VHGT'                  % VECTOR NORMALS %

                                % Offset for entire cell
                                data.land(i_land).hoff = fread(fid,1,'float32');

                                % 65x65 map of (relative!) vertex heights
                                data.land(i_land).vhgt = zeros(65,65,'int8');
                                for ix=1:65
                                for iy=1:65
                                	data.land(i_land).vhgt(ix,iy) = fread(fid,1,'int8');
                                end
                                end
                                clear ix iy;

                                % Junk
                                nul = fread(fid,3,'ubit8');

                                % Continue
                                mem_on = mem_on + 4232;

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            case 'WNAM'                         % MINIMAP %

                            	% World minimap
                                data.land(i_land).wnam = zeros(9,9,'uint8');
                                for ix=1:9
                                for iy=1:9
                                	data.land(i_land).wnam(ix,iy) = fread(fid,1,'uint8');
                                end
                                end
                                clear ix iy;

                                % Continue
                                mem_on = mem_on + 81;

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            case 'VCLR'                   % VERTEX COLORS %

                            	% 65x65 map of RBG values
                                data.land(i_land).vclr = zeros(65,65,3,'uint8');
                                for ix=1:65
                                for iy=1:65
                                	data.land(i_land).vclr(ix,iy,:) = fread(fid,3,'ubit8');
                                end
                                end
                                clear ix iy;

                                % Continue
                                mem_on = mem_on + 12675;

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            case 'VTEX'                 % VERTEX TEXTURES %

                            	% 16x16 map of texture IDs
                                data.land(i_land).vtex = zeros(16,16,'uint16');
                                for ix=1:16
                                for iy=1:16
                                	data.land(i_land).vtex(ix,iy) = fread(fid,1,'uint16');
                                end
                                end
                                clear ix iy;

                                % Continue 
                                mem_on = mem_on + 512;

                            % % % % % % % % % % % % % % % % % % % % % % % %
                            otherwise

                            	%disp([' Skipping LAND sub-entry ' sub_name ' (size ' num2str(sub_size) ')...']);
                                error(['  Unknown LAND sub-entry ' sub_name ' (size ' num2str(sub_size) ')...']);
                                %if(~isempty(sub_size))
                                % 	nul = fread(fid,sub_size,'ubit8');
                                %    mem_on = mem_on + sub_size;
                                %end
                                
                        end %switch sub_name
                        clear sub_*;
                        
                    end %while mem_on < entry_size 
                	clear mem_on;

                %.........................................................%
                otherwise

                	% Skip others for now
                    disp([' Skipping entry ' entry_type ' (size ' num2str(entry_size) ')...']);
                    if(~isempty(entry_size))
                    	nul = fread(fid,entry_size,'ubit8');
                    end
                
            end
            %=============================================================%
            
            % Clean-up
            clear entry_size entry_flags nul;
        
        end %if(~ismember(native2unicode(entry_type),valid_entry_types))
        
        % Next entry
        entry_type = fread(fid,4,'ubit8')';   
        
    end %while(~feof(fid))
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CLEAN DATA OF BAD/EMPTY ENTRIES
    data.file = input_file;
    data = tes3matlab_clean(data);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE DATA TO FILE

    % Save data 
    save(output_file,'data');
    
    % IN FUTURE, WILL BE A NETCDF FOR MORE UNIVERSAL ACCESS
%     % Open file for writing
%     ncid = netcdf.open(output_file,'NETCDF3');
%     
%     % Save header information in file attributes
%     netcdf.putAtt(ncid,0,'file',input_file);
%     netcdf.putAtt(ncid,0,'version',data.hedr.vers);
%     netcdf.putAtt(ncid,0,'type',data.hedr.type);
%     netcdf.putAtt(ncid,0,'company',data.hedr.comp);
%     netcdf.putAtt(ncid,0,'description',data.hedr.desc);
%     
%     % Define dimension lengthss
%     dimid_nr = netcdf.defDim(ncid,'n_regn',i_regn);
%     dimid_nc = netcdf.defDim(ncid,'n_cell',i_cell);
%     dimid_nt = netcdf.defDim(ncid,'n_ltex',i_ltex);
%     dimid_nl = netcdf.defDim(ncid,'n_land',i_land);
%     
%     % Define variables
%     netcdf.defDim(ncid,'regn_
        
end