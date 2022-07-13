function data = tes3matlab_extractdata(input_file,opts,cnst)
% tes3matlab_extractdata(input_file,output_file,opts,cnst)
% Extracts exterior landmass data from TES3 master (.esm) or plugin (.esp)
% file in_file and saves it to the NetCDF file out_file.
% 
% A typical TES3 file contains a list of records, each record containing a
% list of fields.  The first record is always a header (HEDR) record.  The
% HEDR record contains a field representing the number of (remaining) 
% records in the file.
% 
% Mods typically contain MAST records, which come immediately after the
% HEDR but are not included within the HEDR, yet are not included in the
% record count.
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
    % DEFINITIONS
    
    % Variable definitions
    % NOTE: This step is not necessary, but it acts as a sort of glossary
    % for this code.  It may also help when converting this code into other
    % languages which define variables in a preamble (e.g. C languages).
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIT

    % Possible record types
    % NOTE: MAST records are not sub-records of the header, but neither are
    % they counted as a file record, so they must be handled differently...
    rec_types = tes3matlab_recdef;
    %n_types   = numel(rec_types);

    % If file_info is not given, compute it here?
    %file_info = tes3matlab_getfileinfo(input_file);

	% Open the file for reading
	fid = fopen(input_file,'rb');
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HEADER RECORD


	% File signature, also used as name of the Header record
    % (verified it is TES3 in a previous step!)
	header_name = native2unicode(fread(fid,4,'ubit8')');
    if(strcmp(header_name,'TES3')~=true);
        error(['Expected ' input_file ' file header to be TES3, ' ...
               'instead it was ' header_name '!']);
    end
    
    % Get (remainder of) header record size
	header_size = fread(fid,1,'uint32')+8;  % Since there's 8 empty bits
    header_data = fread(fid,header_size,'ubit8');
    file_hedr   = tes3matlab_parserec('TES3',header_data);
    
    % Check that this equals header information if file_info was provided
    % (instead of computed earlier in this function)?
    
    % ...

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OTHER RECORDS
    
    % Init data
    cfll = struct;  i_cell = 0;
    land = struct;  i_land = 0;
    ltex = struct;  i_ltex = 0;
    regn = struct;  i_regn = 0;
    
    % First entry
    i_rec = 1;
    this_type = fread(fid,4,'ubit8')';
    while(~feof(fid))
        
        % Continue if entry_type is valid
        % If it's a MAST entry, its format is a little different from the
        % others, so parse it separately.
        if(~ismember(native2unicode(this_type),rec_types))
            warning([' Unknown entry_type: ' native2unicode(this_type) ', scanning for valid entry_type...']);
        else    
            
            % Record info
            rec_type(i_rec,:) = native2unicode(this_type);
            rec_size(i_rec)   = fread(fid,1,'uint32');
            nul               = fread(fid,4,'ubit8');   % Unused?
            rec_flag(i_rec,:) = fread(fid,4,'ubit8');
            this_data         = fread(fid,rec_size(i_rec),'ubit8');
            this_data         = uint8(this_data);
            
            %=============================================================%
            % Process base on rec_type
            switch native2unicode(this_type)
                
                %---------------------------------------------------------%
                case 'REGN'
                    i_regn = i_regn + 1;
                    [regn(i_regn).name, ...
                     regn(i_regn).fnam, ...
                     regn(i_regn).cnam, ...
                     regn(i_regn).snam, ...
                     regn(i_regn).weat, ...
                     regn(i_regn).bnam  ...
                    ] = tes3matlab_parserec_regn(this_data,cnst.n_weat);
                
                %---------------------------------------------------------%
                case 'LTEX'
                    i_ltex = i_ltex + 1;
                    [ltex(i_ltex).name, ltex(i_ltex).intv, ...
                     ltex(i_ltex).data, ltex(i_ltex).dele ...
                    ] = tes3matlab_parserec_ltex(this_data);
                
                %---------------------------------------------------------%
                case 'CELL'
                    i_cell = i_cell + 1;
                    [cfll(i_cell).name, cfll(i_cell).rgnn, ...
                     cfll(i_cell).flags, cfll(i_cell).x, cfll(i_cell).y,...
                     cfll(i_cell).nam5 ...
                    ] = tes3matlab_parserec_cell(this_data);
                
                %---------------------------------------------------------%
                case 'LAND'
                    i_land = i_land + 1;
                    [land(i_land).dele, land(i_land).x, land(i_land).y, ...
                     land(i_land).flag, land(i_land).vnml, ...
                     land(i_land).hoff, land(i_land).vhgt, ...
                     land(i_land).wnam, land(i_land).vclr, ...
                     land(i_land).vtex ...
                    ] = tes3matlab_parserec_land(this_data,...
                                                 cnst.cell_nh+1,...
                                                 cnst.cell_nt,...
                                                 cnst.cell_nm);
                    % NOTE: # vertices = # cell height faces + 1
                
                %---------------------------------------------------------%
                otherwise;  %warning(['Skipping entry...']);
            
            end
            %==============================================================

            % Next record
            i_rec = i_rec + 1;
            
        end %if(~ismember(native2unicode(this_type),rec_types))
        
        % Next entry
        this_type = fread(fid,4,'ubit8')';   
        
    end %while(~feof(fid))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Debug look?
%     x = [land.x];   xr = [min(x) max(x)];   x0 = xr(1);     nx = diff(xr)+1;
%     y = [land.y];   yr = [min(y) max(y)];   y0 = yr(1);     ny = diff(yr)+1;
%     mapped_wnam = zeros(ny*9,nx*9,'uint8');
%     for i=1:numel(land)
%         ii = double(land(i).x-x0)*9 + (1:9);
%         jj = double(land(i).y-y0)*9 + (1:9);
%         mapped_wnam(jj,ii) = land(i).wnam;
%     end
%     imshow(mapped_wnam(end:-1:1,:));
%     mapped_vtex = zeros(ny*16,nx*16,'uint16');
%     for i=1:numel(land)
%         ii = double(land(i).x-x0)*16 + (1:16);
%         jj = double(land(i).y-y0)*16 + (1:16);
%         mapped_vtex(jj,ii) = land(i).vtex;
%     end
%     imshow(mapped_vtex(end:-1:1,:),[]);
%     mapped_vclr = zeros(ny*65,nx*65,3,'uint8');
%     for i=1:numel(land)
%         ii = double(land(i).x-x0)*65 + (1:65);
%         jj = double(land(i).y-y0)*65 + (1:65);
%         mapped_vclr(jj,ii,:) = land(i).vclr;
%     end
%     imshow(mapped_vclr(end:-1:1,:,:));
%     mapped_ahgt = zeros(ny*65,nx*65);
%     for i=1:numel(land)
%         ii = double(land(i).x-x0)*65 + (1:65);
%         jj = double(land(i).y-y0)*65 + (1:65);
%         mapped_ahgt(jj,ii,:) = tes3matlab_ahgt(land(i).vhgt,land(i).hoff);
%     end
%     imshow(mapped_ahgt(end:-1:1,:),[]);
    
    % Generate and return data structure
    data = struct;
    data.file = input_file;
    data.cell = cfll;
    data.land = land;
    data.ltex = ltex;
    data.regn = regn;
    
end