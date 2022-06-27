function merged_data = tes3matlab_mergedata(opts)

    % Merged data file
    output_file = [opts.dir_head '\data.mat'];
    if(exist(output_file,'file')==2)
        load(output_file);
    else
        
        %=================================================================%
        % INIT
        
        % Get list of data files
        file_list = dir([opts.dir_data '\*.mat']);
        if(isempty(file_list)); error(['No data found in ' opts.dir_data ', run tes3matlab_gatherdata first!']); end

        % Gather datums
        n_file = numel(file_list);
        data   = struct;
        for i=1:n_file
            tmp = load([opts.dir_data '\' file_list(i).name],'data');
            flds = fieldnames(tmp.data);
            if(~isempty(flds))
            for j=1:numel(flds)
                eval(['data(i).' flds{j} '=tmp.data.' flds{j} ';']);
            end
            end
            clear tmp flds;
        end
        clear i;
        
        % Number of records of each data type in each file
        n_regn = zeros(n_file,1,'uint32');
        n_cell = zeros(n_file,1,'uint32');
        n_ltex = zeros(n_file,1,'uint32');
        n_land = zeros(n_file,1,'uint32');
        for i=1:n_file
            if(~isempty(fieldnames(data(i).regn)));    n_regn(i) = numel(data(i).regn);   end
            if(~isempty(fieldnames(data(i).cell)));    n_cell(i) = numel(data(i).cell);   end
            if(~isempty(fieldnames(data(i).ltex)));    n_ltex(i) = numel(data(i).ltex);   end
            if(~isempty(fieldnames(data(i).land)));    n_land(i) = numel(data(i).land);   end
        end
        clear i;
        
        % Init data structures
        merged_data      = struct;
%         flds = fieldnames(data);
%         for i=1:numel(flds)
%             eval(['merged_data.' flds{i} ' = struct;']);
%         end  

        %=================================================================%
        % MERGE REGIONS
        if(any(strcmp(fieldnames(data),'regn')))
            
            % Merge region info
            merged_data.regn = struct;
            i_regn = 1;
            for i=1:n_file
            for j=1:n_regn(i)

                % Region info's
                merged_data.regn(i_regn).name = data(i).regn(j).name;
                merged_data.regn(i_regn).fnam = data(i).regn(j).fnam;
                merged_data.regn(i_regn).cnam = data(i).regn(j).cnam(1:3); % Alpha is always 0
                merged_data.regn(i_regn).weat = data(i).regn(j).weat;
                %merged_data.regn(i_regn).snam = data(i).regn(j).snam;
                %merged_Data.regn(i_regn).bnam = data(i).regn(j).bnam;

                % Older (pre-Bloodmoon) files only have 8 weather entries.
                % Standardize to the modern 10
                if(size(merged_data.regn(i_regn).weat,1)==8)
                    merged_data.regn(i_regn).weat(end+1:end+2) = [0 0];
                end

                % Next region
                i_regn = i_regn+1;

            end
            end
            clear i j i_regn;

            % Determine unique regions
            % NOTE: NAME is what is called by scripts, etc., so go by this
            % NOTE: To follow load order, keep the last of each unique entry so
            %       they overwrite previous entries
            [~,ii]           = unique({merged_data.regn.name},'last');
            merged_data.regn = merged_data.regn(ii);
            nn_regn          = numel(merged_data.regn);

        
        end
        %=================================================================%
        % MERGE CELLS
        if(any(strcmp(fieldnames(data),'cell')))
        
            % Merge cell info
            merged_data.cell = struct;
            i_cell = 1;
            for i=1:n_file
            for j=1:n_cell(i)

                % Region info's
                merged_data.cell(i_cell).name = data(i).cell(j).name;
                merged_data.cell(i_cell).rgnn = data(i).cell(j).rgnn;
               %merged_data.cell(i_cell).data = data(i).cell(j).data;
                merged_data.cell(i_cell).xy   = data(i).cell(j).xy;

                % Next region
                i_cell = i_cell+1;

            end
            end
            clear i j i_cell;

            % Determine unique cells
            % NOTE: Since only exterior cells are saved, can use position
            % NOTE: To follow load order, ensure 'unique' keeps the last value
            [~,ii]           = unique([merged_data.cell.xy]','last','rows');
            merged_data.cell = merged_data.cell(ii);
            nn_cell          = numel(merged_data.cell);

            % Link cells to regions
            if(any(strcmp(fieldnames(data),'regn')))
                regn_names = {merged_data.regn.name};
                for i=1:nn_cell

                    % See if cell's region exists in regn entries
                    ii = find(strcmp(deblank(merged_data.cell(i).rgnn),regn_names)==1);
                    if(~isempty(ii))
                        merged_data.cell(i).regn = ii;
                    end
                    clear ii;

                end
                clear i;
            end
        
        end
        %=================================================================%
        % MERGE LAND TEXTURES
        if(any(strcmp(fieldnames(data),'ltex')))
            
            % Merge ltex info
            i_ltex = 1;
            for i=1:n_file
            for j=1:n_ltex(i)

                % LTEX info's
                merged_data.ltex(i_ltex).file = i;
                merged_data.ltex(i_ltex).name = data(i).ltex(j).name;
                merged_data.ltex(i_ltex).intv = data(i).ltex(j).intv;
                merged_data.ltex(i_ltex).data = data(i).ltex(j).data;

                % Clean 'data' by removing the directory in front of file names
                i1 = [strfind(merged_data.ltex(i_ltex).data,'/'),strfind(merged_data.ltex(i_ltex).data,'\')];	
                if(isempty(i1))	
                    i1 = 0;                                   
                else
                    i1 = max(i1);
                end

                % Clean 'data' by removing the extension
                i2 = strfind(merged_data.ltex(i_ltex).data,'.');         
                if(isempty(i2))
                    i2 = numel(merged_data.ltex(i_ltex).data)+1;
                else
                    i2 = i2(end);
                end

                % Truncate as specified, also clean by setting to lowercase
                ii = (i1+1) : (i2-1);
                merged_data.ltex(i_ltex).data = lower(merged_data.ltex(i_ltex).data(ii));

                % Clean-up
                clear ii i1 i2;

                % Next region
                i_ltex = i_ltex+1;

            end
            end
            clear i j i_ltex;

            % Unique file names
            ltex_data      = {merged_data.ltex.data};
            [ltex_data,aa] = unique(ltex_data);
            nn_ltex        = numel(ltex_data);

            % Don't truncate merged_data.ltex yet, link the intv's first!
            ltex_file = [merged_data.ltex.file];
            intv_link = cell(n_file,0);
            for i=1:n_file

               % This file's data
               ii           = find(ltex_file==i);
               this_intv    = [merged_data.ltex(ii).intv];
               this_data    = {merged_data.ltex(ii).data};
               intv_link{i} = zeros(numel(ii),2);
               clear ii;

               % Link from merged ltex intv to file's original intv
               for j=1:nn_ltex
                   ii = find(strcmp(ltex_data{j},this_data)==1);
                   for k=1:numel(ii)
                       intv_link{i}(ii(k),1) = j;
                       intv_link{i}(ii(k),2) = this_intv(ii(k));
                   end
               end
               clear j ii k;

               % Clean-up
               clear i1 i2 this_* j jj;

            end
            clear i;

            % Save results
            merged_data.intv_link = intv_link;
            merged_data.ltex = merged_data.ltex(aa);
            for i=1:nn_ltex
                merged_data.ltex(i).data = ltex_data{i};
                merged_data.ltex(i).intv = i;
            end
            clear i aa;
        
        end
        %=================================================================%
        % MERGE LANDSCAPES
        if(any(strcmp(fieldnames(data),'land')))
        
            % Merge land info
            merged_data.land = struct;
            i_land = 1;
            for i=1:n_file
            for j=1:n_land(i)
            if(~isempty(data(i).land(j).xy))    

                % Land info's
                merged_data.land(i_land).file = i;
                merged_data.land(i_land).hoff = data(i).land(j).hoff;
                merged_data.land(i_land).vhgt = data(i).land(j).vhgt;
                merged_data.land(i_land).wnam = data(i).land(j).wnam;
                merged_data.land(i_land).xy   = data(i).land(j).xy;

                % Optional entries
                if(isfield(data(i).land(j),'vclr'))
                    merged_data.land(i_land).vclr = data(i).land(j).vclr;
                end
                if(isfield(data(i).land(j),'vnml'))
                    merged_data.land(i_land).vnml = data(i).land(j).vnml;
                end 
                if(isfield(data(i).land(j),'vtex'))
                    merged_data.land(i_land).vtex = data(i).land(j).vtex;
                end

                % Next entry
                i_land = i_land + 1;

            end
            end
            end
            clear i j i_land;

            % Gather the unique land entries
            land_xy          = [merged_data.land.xy]';
            [~,ii]           = unique(land_xy,'rows','last');
            merged_data.land = merged_data.land(ii);
            nn_land          = numel(ii);
            clear ii land_xy; 

            % Process entries
            if(any(strcmp(fieldnames(data),'cell')))
                cell_xy = [merged_data.cell.xy];
                for i=1:nn_land
                    i_cell = find(cell_xy(1,:)==merged_data.land(i).xy(1) & cell_xy(2,:)==merged_data.land(i).xy(2));
                    if(~isempty(i_cell))
                        merged_data.land(i).cell = i_cell;
                        merged_data.land(i).regn = merged_data.cell(i_cell).regn;
                    end
                    clear i_cell;
                end
            end
            if(any(strcmp(fieldnames(data),'ltex')))
            for i=1:nn_land
                
                % Translate vtex to merged_ltex index (if vtex even exists
                if(any(strcmp(fieldnames(data),'ltex')))
                if(isempty(merged_data.land(i).vtex))
                    merged_data.land(i).vtex = zeros(16,16);
                    merged_data.land(i).ltex = zeros(16,16);
                else
                    merged_data.land(i).ltex = zeros(16,16);
                    for ix=1:16
                    for iy=1:16
                    if(merged_data.land(i).vtex(ix,iy)>0)    

                        % NOTE: VTEX LINKS TO INTV+1, SO VTEX=0 MEANS USE _LAND_DEFAULT
                        jj = find(merged_data.intv_link{merged_data.land(i).file}(:,2)+1==merged_data.land(i).vtex(ix,iy));
                        if(numel(jj)==1)
                            merged_data.land(i).ltex(ix,iy) = merged_data.intv_link{merged_data.land(i).file}(jj,1);
                        end   
                    end
                    end
                    end
                    clear ix iy;
                end
                end

                % Translate hoff and vhgt to absolute vertex height (ahgt)
                merged_data.land(i).ahgt = zeros(65,65,'int32');
                if(~isempty(merged_data.land(i).vhgt))
                for ix=1:65
                    merged_data.land(i).ahgt(ix,1) = int32(merged_data.land(i).hoff) + sum(int32(merged_data.land(i).vhgt(1:ix,1)));
                    for iy=2:65
                        merged_data.land(i).ahgt(ix,iy) = merged_data.land(i).ahgt(ix,iy-1) + int32(merged_data.land(i).vhgt(ix,iy));
                    end
                end
                end
                clear ix iy;

            end
            end
            clear i;
        
        end
        %=================================================================%
        % SAVE/BACKUP MERGED DATA
        
        save(output_file,'merged_data');
        
    end
        
end