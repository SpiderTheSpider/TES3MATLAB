function mapped_data = tes3matlab_mapdata(file_data,merged_data,opts,cnst)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIT
    
    % Get size of worldspace
    x = [merged_data.land.x];   xr = [min(x) max(x)];   x0 = xr(1);     nx = diff(xr)+1;
    y = [merged_data.land.y];   yr = [min(y) max(y)];   y0 = yr(1);     ny = diff(yr)+1;

    % Init output data structure
    mapped_data       = struct;
    mapped_data.xr    = xr;
    mapped_data.yr    = yr;
    mapped_data.avail = false(ny,nx);                   % Data availabiliity
    mapped_data.regn  = zeros(ny,nx);                   % Cell region ID
    mapped_data.name  = cell(ny,nx);                    % Cell name
    mapped_data.vtex  = zeros(ny*16,nx*16,'uint16');    % Land texture ID
    mapped_data.vhgt  = -256.*ones(ny*65,nx*65,'single');   % Vertex heights
    mapped_data.vnml  = zeros(ny*65,nx*65,3,'int8');    % Vertex normals
    mapped_data.vclr  = zeros(ny*65,nx*65,3,'uint8');   % Vertex colors
    mapped_data.fhgt  = -256.*ones(ny*64,nx*64,'single');    % Face heights
    mapped_data.fnml  = zeros(ny*64,nx*64,3,'int8');    % Face normals
    mapped_data.fclr  = zeros(ny*64,nx*64,3,'uint8');   % Face colors
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN LOOP
    
    % Loop through landscapes and apply them to the world map
    n_land = numel(merged_data.land);
    for i=1:n_land

        % This land's location on the world map
        ix = merged_data.land(i).x-x0;
        iy = merged_data.land(i).y-y0;
        
        % Data is available at this location
        mapped_data.avail(iy+1,ix+1) = true;

        %=================================================================%
        % MINI-MAP
        
        % Keep the last mini-map loaded for this cell
        this_wnam = merged_data.land(i).wnam;
        if(size(this_wnam,3)>1);    this_wnam = this_wnam(:,:,end); end
        
        % Add to world map
        ii = ix*9 + int32(1:9);
        jj = iy*9 + int32(1:9);
        mapped_data.wnam(jj,ii) = this_wnam;
        
        %=================================================================%
        % CELL REGION ID & NAME
        
        % Keep the last region definition for this cell
        ii = find([merged_data.cell.x]==merged_data.land(i).x & ...
                  [merged_data.cell.y]==merged_data.land(i).y);
        if(~isempty(ii))
            mapped_data.regn(iy+1,ix+1) = merged_data.cell(ii).regn(end);
            mapped_data.name{iy+1,ix+1} = merged_data.cell(ii).rgnn{end};
            if(~isempty(merged_data.cell(ii).name{end}))
                mapped_data.name{iy+1,ix+1} = [mapped_data.name{iy+1,ix+1} ', ' merged_data.cell(ii).name{end}];
            end
        end
        
        %=================================================================%
        % LAND TEXTURE MAP
                
        % Converting to merged_data ltex index
        this_ltex = merged_data.land(i).vtex(:,:,end);
        this_file = merged_data.land(i).file(end);
        glbl_ltex = zeros(size(this_ltex));
        for ii=1:size(this_ltex,1)
        for jj=1:size(this_ltex,2)
        if(~isempty(fieldnames(file_data(this_file).ltex)))
            kk = find( [file_data(this_file).ltex.intv] == this_ltex(ii,jj)-1 );
            if(isempty(kk))
                kk = 0;     % No definition, use _land_default
            else
                kk = find( strcmp({merged_data.ltex.data},file_data(this_file).ltex(kk).data)==1 );
            end
            glbl_ltex(ii,jj) = kk;
        end
        end
        end
        
        % Save to global map
        ii = ix*16 + int32(1:16);
        jj = iy*16 + int32(1:16);
        mapped_data.vtex(jj,ii) = glbl_ltex;
        
        %=================================================================%
        % LANDSCAPE MAPS
        
        % This landscape
        this_vclr = merged_data.land(i).vclr(:,:,:,end);
        this_vnml = merged_data.land(i).vnml(:,:,:,end);
        this_vhgt = tes3matlab_ahgt(merged_data.land(i).vhgt(:,:,end),merged_data.land(i).hoff(end));
        
        % Some junked cells have all depth=0, set to all depth=-256
        if(all(this_vhgt==0))
            this_vhgt = ones(size(this_vhgt),'single').*-256;
        end
        
        % Vertex values
        ii = ix*65 + int32(1:65);
        jj = iy*65 + int32(1:65);
        mapped_data.vclr(jj,ii,:) = this_vclr;
        mapped_data.vnml(jj,ii,:) = this_vnml;
        mapped_data.vhgt(jj,ii,:) = this_vhgt;
        
        % Face values
        ii = ix*64 + int32(1:64);
        jj = iy*64 + int32(1:64);
        mapped_data.fclr(jj,ii,:) = uint8( ( single(this_vclr(1:end-1,1:end-1,:)) ...
                                           + single(this_vclr(2:end,  1:end-1,:)) ...
                                           + single(this_vclr(1:end-1,2:end  ,:)) ...
                                           + single(this_vclr(2:end,  2:end  ,:)) ...
                                           )./4 ...
                                          );
        mapped_data.fnml(jj,ii,:) =  int8( ( single(this_vnml(1:end-1,1:end-1,:)) ...
                                           + single(this_vnml(2:end,  1:end-1,:)) ...
                                           + single(this_vnml(1:end-1,2:end  ,:)) ...
                                           + single(this_vnml(2:end,  2:end  ,:)) ...
                                           )./4 ...
                                          );
        mapped_data.fhgt(jj,ii,:) =      ( ( single(this_vhgt(1:end-1,1:end-1)) ...
                                           + single(this_vhgt(2:end,  1:end-1)) ...
                                           + single(this_vhgt(1:end-1,2:end  )) ...
                                           + single(this_vhgt(2:end,  2:end  )) ...
                                           )./4 ...
                                          );                              
                                      
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SEAM TESTING
    
    % Temporary maps for testing seams in vertex colors and heights
    test_vclr = zeros(ny*65+2,nx*65+2,3);       test_vclr(2:end-1,2:end-1,:) = double(mapped_data.vclr)./255;
    test_vhgt = -256.*ones(ny*65+2,nx*65+2);    test_vhgt(2:end-1,2:end-1)   = mapped_data.vhgt;
    test_vnml = zeros(ny*65+2,nx*65+2,3);       test_vnml(2:end-1,2:end-1,:) = mapped_data.vnml;
    
    % The seam tests: along the north/south lines of the east/west borders
    seamtest_vclr_x = zeros(ny*64,nx+1,3);  
    seamtest_vhgt_x = zeros(ny*64,nx+1);    
    seamtest_vnml_x = zeros(ny*64,nx+1,3); 
    for ix=1:nx+1
        
        ii = (ix-1)*65 + int32(1:2);
        jj = 1:ny*65; jj(65:65:end) = [];
        seamtest_vclr_x(:,ix,:) = uint8(abs(single(test_vclr(jj,ii(1),:)) - single(test_vclr(jj,ii(2),:))));
        seamtest_vnml_x(:,ix,:) = test_vnml(jj,ii(1),:) - test_vclr(jj,ii(2),:);
        seamtest_vhgt_x(:,ix,:) = test_vhgt(jj,ii(1))   - test_vhgt(jj,ii(2));
        
    end
    
    % The seam tests: along the east/west lines of the north/south borders
    seamtest_vclr_y = zeros(ny+1,nx*64,3);
    seamtest_vhgt_y = zeros(ny+1,nx*64);
    seamtest_vnml_y = zeros(ny+1,nx*64,3);
    for iy=1:ny+1
        
        ii = 1:nx*65;   ii(65:65:end) = [];
        jj = (iy-1)*65 + int32(1:2);
        seamtest_vclr_y(iy,:,:) = uint8(abs(single(test_vclr(jj(1),ii,:)) - single(test_vclr(jj(2),ii,:))));
        seamtest_vnml_y(iy,:,:) = test_vnml(jj(1),ii,:) - test_vnml(jj(2),ii,:);
        seamtest_vhgt_y(iy,:,:) = test_vhgt(jj(1),ii)   - test_vhgt(jj(2),ii);
        
        
    end

    % Save seamtests
    mapped_data.seams.vhgt_x = seamtest_vhgt_x;
    mapped_data.seams.vclr_x = seamtest_vclr_x;
    mapped_data.seams.vnml_x = seamtest_vnml_x;
    mapped_data.seams.vhgt_y = seamtest_vhgt_y;
    mapped_data.seams.vclr_y = seamtest_vclr_y;
    mapped_data.seams.vnml_y = seamtest_vnml_y;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FINAL
    
    % Remove row & col #65's
    %mapped_data.vhgt(65:65:end,:)   = [];
    %mapped_data.vhgt(:,65:65:end)   = [];
    %mapped_data.vclr(65:65:end,:,:) = [];
    %mapped_data.vclr(:,65:65:end,:) = [];
    %mapped_data.vnml(65:65:end,:,:) = [];
    %mapped_data.vnml(:,65:65:end,:) = [];
    
    % Compute land-sea mask
    %mapped_data.vmsk = (mapped_data.vhgt>0);
    %mapped_data.fmsk = (mapped_data.fhgt>0);
   
end