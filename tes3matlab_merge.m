function merged_data = tes3matlab_merge(file_data,opts)

    % Load backup if it exists
    bkup_file = [opts.bkup_dir '\' opts.bkup_name '_merged.mat'];
    if(opts.bkup_data==true & exist(bkup_file,'file')==true)
        tmp = load(bkup_file);
        merged_data = tmp.merged_data;
    else

        %=================================================================%
        % INIT

        % Number of data files
        n_file = numel(file_data);

        % Clean empty entries
        for i=1:numel(file_data)
        if(~isempty(fieldnames(file_data(i).ltex)))   
            j = 1;
            while(j<=numel(file_data(i).ltex))
            if(isempty(file_data(i).ltex(j).data))
                file_data(i).ltex(j) = [];
            else
                j=j+1;
            end
            end
        end
        end
        clear i j;

        % Number of records of each data type in each file
        n_regn = zeros(n_file,1,'uint32');
        n_cell = zeros(n_file,1,'uint32');
        n_ltex = zeros(n_file,1,'uint32');
        n_land = zeros(n_file,1,'uint32');
        for i=1:n_file
            if(~isempty(fieldnames(file_data(i).regn)));    n_regn(i) = numel(file_data(i).regn);   end
            if(~isempty(fieldnames(file_data(i).cell)));    n_cell(i) = numel(file_data(i).cell);   end
            if(~isempty(fieldnames(file_data(i).ltex)));    n_ltex(i) = numel(file_data(i).ltex);   end
            if(~isempty(fieldnames(file_data(i).land)));    n_land(i) = numel(file_data(i).land);   end
        end
        clear i;

        % Merged data structure
        merged_data = struct;

        %=================================================================%
        % MERGE REGIONS

        % Merge region info
        i_regn = 1;
        merged_data.regn = struct;
        for i=1:n_file
        for j=1:n_regn(i)

            % Region info's
            merged_data.regn(i_regn).name = file_data(i).regn(j).name;
            merged_data.regn(i_regn).fnam = file_data(i).regn(j).fnam;
            merged_data.regn(i_regn).cnam = file_data(i).regn(j).cnam(1:3); % Alpha is always 0

    % UNUSED FOR NOW, SO COMMENTED OUT TO SAVE PROCESSING TIME
    %         merged_data.regn(i_regn).snam = file_data(i).regn(j).snam;
    %         merged_data.regn(i_regn).weat = file_data(i).regn(j).weat;
    %         
    %         % Older (pre-Bloodmoon) files only have 8 weather entries.
    %         % Standardize to the modern 10
    %         if(size(merged_data.regn(i_regn).weat,1)==8)
    %             merged_data.regn(i_regn).weat(end+1:end+2) = [0 0];
    %         end

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

        % Link to original entries
        regn_names = {merged_data.regn.name};
        for i=1:nn_regn

            % Linkage
            merged_data.regn(i).link = zeros(n_file,1,'uint16');
            for j=1:n_file
            if(~isempty(fieldnames(file_data(j).regn)))

                % See if region exists in this file
                ii = find( strcmp(regn_names{i},{file_data(j).regn.name})==1 );
                if(~isempty(ii))
                    merged_data.regn(i).link(j) = ii;
                end

            end
            end
            clear j ii;

        end
        clear i;

        %=================================================================%
        % MERGE CELLS

        % Merge cell info
        i_cell = 1;
        merged_data.cell = struct;
        for i=1:n_file
        for j=1:n_cell(i)

            % Region info's
            merged_data.cell(i_cell).name = file_data(i).cell(j).name;
            merged_data.cell(i_cell).rgnn = file_data(i).cell(j).rgnn;
           %merged_data.cell(i_cell).data = file_data(i).cell(j).data;
            merged_data.cell(i_cell).xy   = file_data(i).cell(j).xy;

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

        % Link to original entries
        cell_xy = [merged_data.cell.xy];
        for i=1:nn_cell

            % Linkage
            merged_data.cell(i).link = zeros(n_file,1,'uint16');
            for j=1:n_file
            if(isfield(file_data(j).cell,'xy'))

                % This file's positions
                this_xy = [file_data(j).cell.xy];

                % See if cell exists in this file
                ii = find( cell_xy(1,:)==this_xy(1) & cell_xy(2,:)==this_xy(2) );
                if(~isempty(ii))
                    merged_data.cell(i).link(j) = ii;
                end

            end
            end
            clear j ii this_xy;

        end
        clear i;

        % Link to regions
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

        %=================================================================%
        % MERGE LAND TEXTURES

        % Merge ltex info
        i_ltex = 1;
        merged_data.ltex = struct;
        for i=1:n_file
        for j=1:n_ltex(i)

            % Region info's
            merged_data.ltex(i_ltex).file = i;
            merged_data.ltex(i_ltex).name = file_data(i).ltex(j).name;
            merged_data.ltex(i_ltex).intv = file_data(i).ltex(j).intv;
            merged_data.ltex(i_ltex).data = file_data(i).ltex(j).data;

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

        %=================================================================%
        % MERGE LANDSCAPES

        % Merge land info
        i_land = 1;
        merged_data.land = struct;
        for i=1:n_file
        for j=1:n_land(i)
        if(~isempty(file_data(i).land(j).xy))    

            % Land info's
            merged_data.land(i_land).file = i;
            merged_data.land(i_land).hoff = file_data(i).land(j).hoff;
            %merged_land(i_land).vclr = file_data(i).land(j).vclr;
            merged_data.land(i_land).vhgt = file_data(i).land(j).vhgt;
            %merged_land(i_land).vnml = file_data(i).land(j).vnml;
            merged_data.land(i_land).wnam = file_data(i).land(j).wnam;
            merged_data.land(i_land).xy   = file_data(i).land(j).xy;

            % Sometimes vtex doesn't exist even when vhgt, etc. do
            if(isfield(file_data(i).land(j),'vtex'))
                merged_data.land(i_land).vtex = file_data(i).land(j).vtex;
            else
                merged_data.land(i_land).vtex = [];
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
        cell_xy = [merged_data.cell.xy];
        for i=1:nn_land

            % Link to cell
            i_cell = find(cell_xy(1,:)==merged_data.land(i).xy(1) & cell_xy(2,:)==merged_data.land(i).xy(2));
            if(~isempty(i_cell))
                merged_data.land(i).cell = i_cell;
                merged_data.land(i).regn = merged_data.cell(i_cell).regn;
            end
            clear i_cell;

            % Translate vtex to merged_ltex index (if vtex even exists
            if(isempty(merged_data.land(i).vtex))
                merged_data.land(i).vtex = zeros(16,16);
                merged_data.land(i).ltex = zeros(16,16);
            else
                merged_data.land(i).ltex = zeros(16,16);
                for ix=1:16
                for iy=1:16
                if(merged_data.land(i).vtex(ix,iy)>0)    

                    jj = find(merged_data.intv_link{merged_data.land(i).file}(:,2)+1==merged_data.land(i).vtex(ix,iy));
                    if(numel(jj)~=1)
                        pause(1);
                    else
                        merged_data.land(i).ltex(ix,iy) = merged_data.intv_link{merged_data.land(i).file}(jj,1);
                    end

                end
                end
                end
                clear ix iy;
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

    %         % Compute face height (fhgt) from absolute vertex height
    %         tmp = single(merged_data.land(i).ahgt);
    %         merged_data.land(i).fhgt = (tmp(1:end-1,1:end-1) + tmp(2:end,1:end-1) + tmp(1:end-1,2:end) + tmp(2:end,2:end) )./4;
    %         clear tmp;

        end
        clear i;

        % Backup if requested
        if(opts.bkup_data==true)
            save(bkup_file,'merged_data');
        end
    
    end

end