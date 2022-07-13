function merged_data = tes3matlab_mergedata(file_data,opts,cnst)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIT
    
    % Number of data files
    n_file = numel(file_data);
    
    % Number of records of each data type in each file
    fld_name = {'regn','cell','ltex','land'};
    n_flds   = numel(fld_name);
    n_fld    = zeros(n_file,n_flds);
    for i=1:n_file
    for j=1:n_flds
        eval(['tmp=file_data(i).' fld_name{j} ';']);
        if(~isempty(tmp));  n_fld(i,j) = numel(tmp);  end
    end
    end
    clear i;

    % Init merged data structure
    merged_data = struct;
    merged_data.file = {file_data.file};
    merged_data.regn = struct;
     merged_data.regn.name = '';
    merged_data.cell = struct;
     merged_data.cell.x = [];
     merged_data.cell.y = [];
    merged_data.ltex = struct;
     merged_data.ltex.data = '';
    merged_data.land = struct;
     merged_data.land.x = [];
     merged_data.land.y = [];
    %for i=1:n_flds; eval(['merged_data.' fld_name{i} ' = struct;']);    end
    
    % Region indices
    i_fld = zeros(n_flds,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MERGE REGIONS
    if(any(strcmp( fieldnames(file_data),'regn' )))
        
        % Link to global field info's
        ii = find(strcmp(fld_name,'regn')==1);
        
        % Gather region info
        for i=1:n_file
        for j=1:n_fld(i,ii)  
        if(~isempty(fieldnames(file_data(i).regn(j))))
            
            % Determine if this is a new region by its NAME field
            jj = find(strcmp({merged_data.regn.name},file_data(i).regn(j).name));
            if(isempty(jj))

                % A new region!
                i_fld(ii) = i_fld(ii)+1;
                
                % Region info's
                merged_data.regn(i_fld(ii)).name = file_data(i).regn(j).name;
                merged_data.regn(i_fld(ii)).fnam = {file_data(i).regn(j).fnam};
                merged_data.regn(i_fld(ii)).cnam = file_data(i).regn(j).cnam(1:3);  % 4th value, alpha, is always 0
                merged_data.regn(i_fld(ii)).weat = file_data(i).regn(j).weat;
                merged_data.regn(i_fld(ii)).file = i;
                merged_data.regn(i_fld(ii)).n    = 1;
            
            else
                
                % In case region info's are changed, save as new entries
                merged_data.regn(jj).fnam{end+1}   = file_data(i).regn(j).fnam;
                merged_data.regn(jj).cnam(:,end+1) = file_data(i).regn(j).cnam(1:3);
                merged_data.regn(jj).weat(end+1,:) = file_data(i).regn(j).weat;
                merged_data.regn(jj).file(end+1)   = i;
                merged_data.regn(jj).n             = merged_data.regn(jj).n + 1;

            end
            
        end
        end
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MERGE CELLS
    if(any(strcmp( fieldnames(file_data),'cell')))

        % Link to global field info's
        ii = find(strcmp(fld_name,'cell')==1);
        
        % Gather region info
        for i=1:n_file
        for j=1:n_fld(i,ii)    
        if(~isempty(fieldnames(file_data(i).cell(j))))
        if(~isempty(file_data(i).cell(j).rgnn)) % If is exterior cell...
            
            % Determine if this is a new cell by its coordinates
            jj = find([merged_data.cell.x]==file_data(i).cell(j).x & ...
                      [merged_data.cell.y]==file_data(i).cell(j).y);
            if(isempty(jj))
                
                % A new cell!
                i_fld(ii) = i_fld(ii)+1;
                
                % Save cell info
                merged_data.cell(i_fld(ii)).name  = {file_data(i).cell(j).name};
                merged_data.cell(i_fld(ii)).rgnn  = {file_data(i).cell(j).rgnn};
                merged_data.cell(i_fld(ii)).flags = file_data(i).cell(j).flags;
                merged_data.cell(i_fld(ii)).nam5  = file_data(i).cell(j).nam5(:);
                merged_data.cell(i_fld(ii)).x     = file_data(i).cell(j).x;
                merged_data.cell(i_fld(ii)).y     = file_data(i).cell(j).y;
                merged_data.cell(i_fld(ii)).file  = i;
                merged_data.cell(i_fld(ii)).n     = 1;
                
                % Link to region
                merged_data.cell(i_fld(ii)).regn = find( strcmp({merged_data.regn.name},merged_data.cell(i_fld(ii)).rgnn)==1 );
                
            else
                
                % Append
                merged_data.cell(jj).name{end+1}    = file_data(i).cell(j).name;
                merged_data.cell(jj).rgnn{end+1}    = file_data(i).cell(j).rgnn;
                merged_data.cell(jj).flags(:,end+1) = file_data(i).cell(j).flags;
                merged_data.cell(jj).nam5(:,end+1)  = file_data(i).cell(j).nam5(:);
                merged_data.cell(jj).file(end+1)    = i;
                merged_data.cell(jj).n              = merged_data.cell(jj).n + 1;
                
                % Link to region
                merged_data.cell(jj).regn(end+1) = find( strcmp({merged_data.regn.name},merged_data.cell(jj).rgnn{end})==1 );

            end
        
        end
        end
        end
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MERGE LAND TEXTURES
    if(any(strcmp( fieldnames(file_data),'ltex')))
        
        % Link to global field info's
        ii = find(strcmp(fld_name,'ltex')==1);
        
        % Gather region info
        for i=1:n_file
        for j=1:n_fld(i,ii)    
        if(~isempty(fieldnames(file_data(i).ltex(j))))
        if(~file_data(i).ltex(j).dele)  % Don't merge if a 'deleted' texture...
            
            % Determine if this is a new ltex by its file name
            jj = find(strcmp({merged_data.ltex.data},file_data(i).ltex(j).data)==1);
            if(isempty(jj))
                
                % A new texture!
                i_fld(ii) = i_fld(ii) + 1;
                
                % Save land texture data
                merged_data.ltex(i_fld(ii)).name = {file_data(i).ltex(j).name};
                merged_data.ltex(i_fld(ii)).intv = file_data(i).ltex(j).intv;
                merged_data.ltex(i_fld(ii)).data = file_data(i).ltex(j).data;
                merged_data.ltex(i_fld(ii)).file = i;
                merged_data.ltex(i_fld(ii)).n    = 1;
                
            else
                
                % Append
                merged_data.ltex(jj).name{end+1} = file_data(i).ltex(j).name;
                merged_data.ltex(jj).intv(end+1) = file_data(i).ltex(j).intv;
                merged_data.ltex(jj).file(end+1) = i;
                merged_data.ltex(jj).n           = merged_data.ltex(jj).n + 1;

            end

        end
        end
        end
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MERGE LANDSCAPES
    if(any(strcmp( fieldnames(file_data),'land')))
        
        % Link to global field info's
        ii = find(strcmp(fld_name,'land')==1);
        
        % Gather region info
        for i=1:n_file
        if(opts.merge_land(i))
        for j=1:n_fld(i,ii)    
        if(~isempty(fieldnames(file_data(i).land(j))))
        if(~file_data(i).land(j).dele)
            
            % Determine if this is a new cell by its coordinates
            jj = find([merged_data.land.x]==file_data(i).land(j).x & ...
                      [merged_data.land.y]==file_data(i).land(j).y);
            if(isempty(jj))
                
                % A new landscape!
                i_fld(ii) = i_fld(ii)+1;
                
                % Save land info
                merged_data.land(i_fld(ii)).x    = file_data(i).land(j).x;
                merged_data.land(i_fld(ii)).y    = file_data(i).land(j).y;
                merged_data.land(i_fld(ii)).flag = file_data(i).land(j).flag;
                merged_data.land(i_fld(ii)).vnml = file_data(i).land(j).vnml;
                merged_data.land(i_fld(ii)).hoff = file_data(i).land(j).hoff;
                merged_data.land(i_fld(ii)).vhgt = file_data(i).land(j).vhgt;
                merged_data.land(i_fld(ii)).wnam = file_data(i).land(j).wnam;
                merged_data.land(i_fld(ii)).vclr = file_data(i).land(j).vclr;
                merged_data.land(i_fld(ii)).vtex = file_data(i).land(j).vtex;
                merged_data.land(i_fld(ii)).file = i;
                merged_data.land(i_fld(ii)).n    = 1;
                
            else
                
                % Append land info
                merged_data.land(jj).flag(:,end+1)     = file_data(i).land(j).flag;
                merged_data.land(jj).vnml(:,:,:,end+1) = file_data(i).land(j).vnml;
                merged_data.land(jj).hoff(end+1)       = file_data(i).land(j).hoff;
                merged_data.land(jj).vhgt(:,:,end+1)   = file_data(i).land(j).vhgt;
                merged_data.land(jj).wnam(:,:,end+1)   = file_data(i).land(j).wnam;
                merged_data.land(jj).vclr(:,:,:,end+1) = file_data(i).land(j).vclr;
                merged_data.land(jj).vtex(:,:,end+1)   = file_data(i).land(j).vtex;
                merged_data.land(jj).file(end+1)       = i; 
                merged_data.land(jj).n                 = merged_data.land(jj).n + 1;
                
            end
            
        end
        end
        end
        end
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
end