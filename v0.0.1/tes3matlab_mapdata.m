function map_data = tes3matlab_mapdata(merged_data,opts)

    % If map_list = {'ALL'}, set all maps
    if(numel(opts.map_list)==1 && strcmp(opts.map_list{1},'ALL'))
        opts.map_list = {'MASK','WNAM','REGN','CELL','LTEX','VHGT','FHGT','VNML','VCLR','BLOK'};
    end

    % Dimensions
    n_land = numel(merged_data.land);
    n_cell = numel(merged_data.cell);

    % Get world coordinates and limits
    world_xy = [merged_data.land.xy];
    xr = [min(world_xy(1,:)), max(world_xy(1,:))];  x0 = xr(1); nx = diff(xr)+1;
    yr = [min(world_xy(2,:)), max(world_xy(2,:))];  y0 = yr(1); ny = diff(yr)+1;
    
    % Init output data structure and save world limits
    map_data = struct;
    map_data.x = xr;
    map_data.y = yr;

    % Loop through maps to generate
    for z=1:numel(opts.map_list)
    switch opts.map_list{z}

        %=================================================================%
        % The MASK map identifies which datums are available
        case 'MASK'
            
            % Masks
            map_data.mask = false(nx,ny,5);
            for i=1:n_land
                
                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;

                % Has wnam data
                if(~isempty(merged_data.land(i).wnam))
                if(~all(merged_data.land(i).wnam(:) == merged_data.land(i).wnam(1)))
                    map_data.mask(ix,iy,1) = 1;
                end
                end
                
                % Has cell and region data
                if(~isempty(merged_data.land(i).cell))
                if(~isempty(merged_data.cell(merged_data.land(i).cell).rgnn))
                    map_data.mask(ix,iy,2) = 1;
                end
                if(~isempty(merged_data.cell(merged_data.land(i).cell).regn))
                    map_data.mask(ix,iy,3) = 1;
                end
                end
                
                % has land data, has land data above z=0
                if(~isempty(merged_data.land(i).ahgt))
                if(~all(merged_data.land(i).ahgt(:)==merged_data.land(i).ahgt(1)))
                	map_data.mask(ix,iy,4) = 1;
                    if(any(merged_data.land(i).ahgt(:)>0))
                        map_data.mask(ix,iy,5) = 1;
                    end
                end
                end
                
                % If land textures are available
                if(~isempty(merged_data.land(i).vtex))
                    map_data.mask(ix,iy,6) = 1;
                end

                % Clean-up
                clear ix iy ii jj tmp;
                
            end
            map_data.mask = permute(map_data.mask,[2 1 3]);
            clear i;
            
        %=================================================================%
        % WORLD MINIMAP
        case 'WNAM'
    
            % The world mini-map
            map_data.wnam = zeros(nx*9,ny*9);
            for i=1:n_land
            if(~isempty(merged_data.land(i).wnam))
                
                % Location
                ix = merged_data.land(i).xy(1)-x0+1;   
                iy = merged_data.land(i).xy(2)-y0+1; 

                % Subspace
                ii = (ix-1)*9+1 : ix*9;
                jj = (iy-1)*9+1 : iy*9;
                map_data.wnam(ii,jj) = merged_data.land(i).wnam'; 

            end
            end
            map_data.wnam = map_data.wnam';
            clear i ix iy ii jj;
            
        %=================================================================%
        % REGION MAP
        case 'REGN'
            
            % The world region map
            map_data.regn = zeros(nx,ny);
            map_data.rclr = zeros(nx,ny,3,'uint8');
            for i=1:n_cell
            if(~isempty(merged_data.cell(i).regn))
            %if(all(abs(merged_data.cell(i).xy)<255))
                
                % Location
                ix = merged_data.cell(i).xy(1)-x0+1;
                iy = merged_data.cell(i).xy(2)-y0+1;
                
                % Sometimes interior cells or bugged cells have HUGE xy's
                if((ix+x0-1)>=xr(1) & (iy+y0-1)>=yr(1) & ...
                   (ix+x0-1)<=xr(2) & (iy+y0-1)<=yr(2) )

                    % Region
                    map_data.regn(ix,iy) = merged_data.cell(i).regn;
                    map_data.rclr(ix,iy,:) = merged_data.regn(merged_data.cell(i).regn).cnam;
                    
                end
                clear ix iy;

            %end
            end
            end
            map_data.regn = map_data.regn;
            map_data.rclr = permute(map_data.rclr,[2 1 3]);
            clear i ix iy;
            
        %=================================================================%
        % CELL MAP
        case 'CELL'
            
            % The world cell map
            map_data.cell = cell(nx,ny,4);
            for i=1:n_cell

                % Location
                ix = merged_data.cell(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.cell(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Sometimes interior cells or bugged cells have HUGE xy's
                if((ix+x0-1)>=xr(1) & (iy+y0-1)>=yr(1) & ...
                   (ix+x0-1)<=xr(2) & (iy+y0-1)<=yr(2) )
                            
                    % Cell name
                    map_data.cell{ix,iy,1} = merged_data.cell(i).name;                   % Point of interest
                    if(isempty(merged_data.cell(i).regn))
                        map_data.cell{ix,iy,2} = merged_data.cell(i).rgnn;
                        map_data.cell{ix,iy,3} = merged_data.cell(i).rgnn;
                    else
                        map_data.cell{ix,iy,2} = merged_data.regn(merged_data.cell(i).regn).name; % Subregion
                        map_data.cell{ix,iy,3} = merged_data.regn(merged_data.cell(i).regn).fnam; % Region+Subregion
                    end
                    
                end
                clear ix iy;

            end
            map_data.cell = permute(map_data.cell,[2 1 3]);
            clear i ix iy;
            
        %=================================================================%
        % TEXTURE MAP
        case 'LTEX'
            
            % The world texture map
            map_data.ltex = zeros(nx*16,ny*16);
            for i=1:n_land

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Subspace
                ii = (ix-1)*16+1 : ix*16;
                jj = (iy-1)*16+1 : iy*16;

                % LTEX
                map_data.ltex(ii,jj) = tes3matlab_translate_vtex(merged_data.land(i).ltex);

            end
            map_data.ltex = map_data.ltex';
            clear i ix iy ii jj;
            
        %=================================================================%
        % VERTEX HEIGHT MAP
        case 'VHGT'  
            
            % The world vertex height map
            map_data.vhgt = ones(nx*65,ny*65).*-256;    % Vertex heights
            for i=1:n_land

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Subspace
                ii = (ix-1)*65+1 : ix*65;
                jj = (iy-1)*65+1 : iy*65;

                % HGHT
                map_data.vhgt(ii,jj) = merged_data.land(i).ahgt';             

            end
            map_data.vhgt = map_data.vhgt';
            clear i ix iy ii jj;
            
        %=================================================================%
        % VERTEX COLOR MAP
        case 'VCLR'  
            
            % The world vertex height map
            map_data.vclr = ones(nx*65,ny*65,3,'uint8');    % Vertex heights
            for i=1:n_land
            if(~isempty(merged_data.land(i).vclr))    

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Subspace
                ii = (ix-1)*65+1 : ix*65;
                jj = (iy-1)*65+1 : iy*65;

                % HGHT
                map_data.vclr(ii,jj,:) = permute(merged_data.land(i).vclr,[2 1 3]);             

            end
            end
            map_data.vclr = permute(map_data.vclr,[2 1 3]);
            clear i ix iy ii jj;
            
        %=================================================================%
        % VERTEX NORMAL MAP
        case 'VNML'  
            
            % The world vertex normals map
            map_data.vnml = ones(nx*65,ny*65,3);    % Vertex heights
            for i=1:n_land
            if(~isempty(merged_data.land(i).vnml))

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;
                iy = merged_data.land(i).xy(2)-y0+1;
                
                % Subspace
                ii = (ix-1)*65+1 : ix*65;
                jj = (iy-1)*65+1 : iy*65;

                % HGHT
                map_data.vnml(ii,jj,:) = permute(merged_data.land(i).vnml,[2 1 3]);             

            end
            end
            map_data.vnml = permute(map_data.vnml,[2 1 3]);
            clear i ix iy ii jj;
            
        %=================================================================%
        % FACE HEIGHT MAP
        case 'FHGT'
            
            % The world cell map
            map_data.fhgt = ones(nx*64,ny*64).*-256;    % Face heights
            for i=1:n_land

                % Location
                ix = merged_data.land(i).xy(1)-x0+1;    %ix = (ix-1)*9+1;
                iy = merged_data.land(i).xy(2)-y0+1;   % iy = (iy-1)*9+1;
                
                % Subspace
                ii = (ix-1)*64+1 : ix*64;
                jj = (iy-1)*64+1 : iy*64;

                % HGHT
                map_data.fhgt(ii,jj) = ( merged_data.land(i).ahgt(1:64,1:64) + ...
                                     merged_data.land(i).ahgt(1:64,2:65) + ...
                                     merged_data.land(i).ahgt(2:65,1:64) + ...
                                     merged_data.land(i).ahgt(2:65,2:65) )./4;
                map_data.fhgt(ii,jj) = map_data.fhgt(ii,jj)';                 

            end
            map_data.fhgt = map_data.fhgt';
            clear i ix iy ii jj;
            
        %=================================================================%
        % BLOCK MAP
        case 'BLOK'
            
            % The world block height map
            nn = opts.dims.cell_ny;
            zz = (opts.dims.cell_ng/opts.dims.hinc_ng)/nn;
            ff = 64/nn;
            
            % Interpolate face height onto new map
            [y1,x1] = meshgrid(0:(nx*64-1),0:(ny*64-1));
            [y2,x2] = meshgrid(0:ff:nx*64-ff,0:ff:ny*64-ff);
            map_data.bhgt = interp2(y1,x1,map_data.fhgt,y2,x2);
            
            % Round to nearest blok unit
            map_data.bhgt = ceil(map_data.bhgt./zz);
            
        %=================================================================%
        otherwise;  warning(['Unknown map type: ' opts.map_list{z} ', skipping...']);
            
    end
    end
    clear z;
    
    % Save?
    save([opts.dir_head '\maps.mat'],'map_data');

end