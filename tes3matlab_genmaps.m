function maps = tes3matlab_genmaps(merged_data,opts)

    % If map_list = {'ALL'}, set all maps
    if(numel(opts.map_list)==1 && strcmp(opts.map_list{1},'ALL'))
        opts.map_list = {'MASK','WNAM','REGN','CELL','LTEX','VHGT','FHGT'};
    end

    % Dimensions
    n_land = numel(merged_data.land);
    n_cell = numel(merged_data.cell);

    % Get world coordinates and limits
    world_xy = [merged_data.land.xy];
    xr = [min(world_xy(1,:)), max(world_xy(1,:))];  x0 = xr(1); nx = diff(xr)+1;
    yr = [min(world_xy(2,:)), max(world_xy(2,:))];  y0 = yr(1); ny = diff(yr)+1;
    
    % Init output data structure and save world limits
    maps = struct;
    maps.x = xr;
    maps.y = yr;

    % Loop through maps to generate
    for z=1:numel(opts.map_list)
    switch opts.map_list{z}
        
        %=================================================================%
        % The MASK map identifies which datums are available
        case 'MASK'
            
            % Masks
            maps.mask = false(nx,ny,5);
            for i=1:n_land
                
                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;

                % Has wnam data
                if(~isempty(merged_data.land(i).wnam))
                if(~all(merged_data.land(i).wnam(:) == merged_data.land(i).wnam(1)))
                    maps.mask(ix,iy,1) = 1;
                end
                end
                
                % Has cell and region data
                if(~isempty(merged_data.land(i).cell))
                if(~isempty(merged_data.cell(merged_data.land(i).cell).rgnn))
                    maps.mask(ix,iy,2) = 1;
                end
                if(~isempty(merged_data.cell(merged_data.land(i).cell).regn))
                    maps.mask(ix,iy,3) = 1;
                end
                end
                
                % Percentage of points above hght==0
                if(~isempty(merged_data.land(i).ahgt))
                if(~all(merged_data.land(i).ahgt(:)==merged_data.land(i).ahgt(1)))
                	maps.mask(ix,iy,4) = 1;
                    if(any(merged_data.land(i).ahgt(:)>0))
                        maps.mask(ix,iy,5) = 1;
                    end
                end
                end
                
                % If land textures are available
                if(~isempty(merged_data.land(i).vtex))
                    maps.mask(ix,iy,6) = 1;
                end

                % Clean-up
                clear ix iy ii jj tmp;
                
            end
            maps.mask = permute(maps.mask,[2 1 3]);
            clear i;

        %=================================================================%
        % WORLD MINIMAP
        case 'WNAM'
    
            % The world mini-map
            maps.wnam = zeros(nx*9,ny*9);
            for i=1:n_land
            if(~isempty(merged_data.land(i).wnam))
                
                % Location
                ix = merged_data.land(i).xy(1)-x0+1;   
                iy = merged_data.land(i).xy(2)-y0+1; 

                % Subspace
                ii = (ix-1)*9+1 : ix*9;
                jj = (iy-1)*9+1 : iy*9;
                maps.wnam(ii,jj) = merged_data.land(i).wnam'; 

            end
            end
            maps.wnam = maps.wnam';
            clear i ix iy ii jj;
       
        %=================================================================%
        % REGION MAP
        case 'REGN'
            
            % The world region map
            maps.regn = zeros(nx,ny,3,'uint8');
            for i=1:n_cell
            if(~isempty(merged_data.cell(i).regn))
            if(all(abs(merged_data.cell(i).xy)<255))
                
                % Location
                ix = merged_data.cell(i).xy(1)-x0+1;
                iy = merged_data.cell(i).xy(2)-y0+1;
                
                % Sometimes interior cells or bugged cells have HUGE xy's
                if((ix+x0-1)>=xr(1) & (iy+y0-1)>=yr(1) & ...
                   (ix+x0-1)<=xr(2) & (iy+y0-1)<=yr(2) )

                    % Region
                    maps.regn(ix,iy,:) = merged_data.regn(merged_data.cell(i).regn).cnam;
                    
                end
                clear ix iy;

            end
            end
            end
            maps.regn = permute(maps.regn,[2 1 3]);
            clear i ix iy;
            
        %=================================================================%
        % CELL MAP
        case 'CELL'
            
            % The world cell map
            maps.cell = cell(nx,ny,4);
            for i=1:n_cell

                % Location
                ix = merged_data.cell(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.cell(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Sometimes interior cells or bugged cells have HUGE xy's
                if((ix+x0-1)>=xr(1) & (iy+y0-1)>=yr(1) & ...
                   (ix+x0-1)<=xr(2) & (iy+y0-1)<=yr(2) )
                            
                    % Cell name
                    maps.cell{ix,iy,1} = merged_data.cell(i).name;                   % Point of interest
                    if(isempty(merged_data.cell(i).regn))
                        maps.cell{ix,iy,2} = merged_data.cell(i).rgnn;
                        maps.cell{ix,iy,3} = merged_data.cell(i).rgnn;
                    else
                        maps.cell{ix,iy,2} = merged_data.regn(merged_data.cell(i).regn).name; % Subregion
                        maps.cell{ix,iy,3} = merged_data.regn(merged_data.cell(i).regn).fnam; % Region+Subregion
                    end
                    
                end
                clear ix iy;

            end
            maps.cell = permute(maps.cell,[2 1 3]);
            clear i ix iy;
            
        %=================================================================%
        % TEXTURE MAP
        case 'LTEX'
            
            % The world texture map
            maps.ltex = zeros(nx*16,ny*16);
            for i=1:n_land

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Subspace
                ii = (ix-1)*16+1 : ix*16;
                jj = (iy-1)*16+1 : iy*16;

                % LTEX
                maps.ltex(ii,jj) = tes3matlab_translate_vtex(merged_data.land(i).ltex);

            end
            maps.ltex = maps.ltex';
            clear i ix iy ii jj;
            
        %=================================================================%
        % VERTEX HEIGHT MAP
        case 'VHGT'  
            
            % The world vertex height map
            maps.vhgt = ones(nx*65,ny*65).*-256;    % Vertex heights
            for i=1:n_land

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Subspace
                ii = (ix-1)*65+1 : ix*65;
                jj = (iy-1)*65+1 : iy*65;

                % HGHT
                maps.vhgt(ii,jj) = merged_data.land(i).ahgt';             

            end
            maps.vhgt = maps.vhgt';
            clear i ix iy ii jj;
            
        %=================================================================%
        % FACE HEIGHT MAP
        case 'FHGT'
            
            % The world cell map
            maps.fhgt = ones(nx*64,ny*64).*-256;    % Face heights
            for i=1:n_land

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Subspace
                ii = (ix-1)*64+1 : ix*64;
                jj = (iy-1)*64+1 : iy*64;

                % HGHT
                maps.fhgt(ii,jj) = ( merged_data.land(i).ahgt(1:64,1:64) + ...
                                     merged_data.land(i).ahgt(1:64,2:65) + ...
                                     merged_data.land(i).ahgt(2:65,1:64) + ...
                                     merged_data.land(i).ahgt(2:65,2:65) )./4;
                maps.fhgt(ii,jj) = maps.fhgt(ii,jj)';                 

            end
            maps.fhgt = maps.fhgt';
            clear i ix iy ii jj;

        %=================================================================%
        otherwise;  warning(['Unknown map type: ' opts.map_list{z} ', skipping...']);
        
    end
    end

end