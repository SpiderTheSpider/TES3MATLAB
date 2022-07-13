function imaged_data = tes3matlab_drawdata(mapped_data,merged_data,opts,cnst)
% 
% Draws raw-scale images of the following data:
% 1: RAW
% > 01	bmap    basemap             064x064
% > 02  mmap    minimap             009x009
% > 03  regn    region_colors       001x001
% > 04  vtex    texture_indices     016x016
% > 05  ltex    texture_colors      016x016
% > 06  vclr    vertex_colors       064x064
% > 07  vnml    vertex_normals      064x064
% > 08  fhgt    face_heights        064x064
% 2: DERIVED
% 

% NOTES
% Currently, world min/max is
%  > [-232.250 330.750] (TR, YARDS)
%  > [-291.125 678.125] (PT/TR, YARDS)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIT
    
%     % World coordinates
%     if(~isempty(opts.draw_xr))
%         x = opts.draw_xr(1) : opts.draw_xr(2);
%     else
%         x = min([file_data.land.x]) : max([file_data.land.x]);
%     end
%     if(~isempty(opts.draw_yr))
%         y = opts.draw_yr(1) : opts.draw_yr(2);
%     else
%         x = min([file_data.land.y]) : max([file_data.land.y]);
%     end
%     
%     % World cell range & size
%     xr = [x(1) x(end)];     nx = x(end)-x(1)+1;     x0 = min([file_data.land.x])-x(1)+1;
%     yr = [y(1) y(end)];     ny = y(end)-y(1)+1;     y0 = min([file_data.land.y])-y(1)+1;
%     
    % World dimensions
    [ny,nx] = size(mapped_data.avail);
    ng      = opts.seam_grid;
    
    % Maps structure
    imaged_data = struct;
    imaged_data.bmap.rgb = zeros(ny*064,nx*064,3,'uint8');
    imaged_data.bmap.alp = zeros(ny*064,nx*064,1,'uint8');
    imaged_data.mmap.rgb = zeros(ny*009,nx*009,3,'uint8');
    imaged_data.mmap.alp = zeros(ny*009,nx*009,1,'uint8');
    imaged_data.regn.rgb = zeros(ny*001,nx*001,3,'uint8');
    imaged_data.regn.alp = zeros(ny*001,nx*001,1,'uint8');
    imaged_data.vtex.rgb = zeros(ny*16,nx*16,3,'uint8');
    imaged_data.vtex.alp = zeros(ny*16,nx*16,1,'uint8');
    imaged_data.lclr.rgb = zeros(ny*16,nx*16,3,'uint8');
    imaged_data.lclr.alp = zeros(ny*16,nx*16,1,'uint8');
    imaged_data.vclr.rgb = zeros(ny*64,nx*64,3,'uint8');
    imaged_data.vclr.alp = zeros(ny*64,nx*64,1,'uint8');
    imaged_data.vnml.rgb = zeros(ny*64,nx*64,3,'uint8');
    imaged_data.vnml.alp = zeros(ny*64,nx*64,1,'uint8');
    imaged_data.fhgt.rgb = zeros(ny*64,nx*64,3,'single');
    imaged_data.fhgt.alp = zeros(ny*64,nx*64,1,'uint8');
    imaged_data.tshd.rgb = zeros(ny*64,nx*64,3,'uint8');
    imaged_data.tshd.alp = zeros(ny*64,nx*64,1,'uint8');
    imaged_data.tcnt.rgb = zeros(ny*64,nx*64,3,'uint8');
    imaged_data.tcnt.alp = zeros(ny*64,nx*64,1,'uint8');
    imaged_data.step.rgb = zeros(ny*64,nx*64,3,'uint8');
    imaged_data.step.alp = zeros(ny*64,nx*64,1,'uint8');
    imaged_data.seam_vhgt.rgb = zeros(ny*ng,nx*ng,3,'uint8');
    imaged_data.seam_vhgt.alp = zeros(ny*ng,nx*ng,1,'uint8');
    imaged_data.seam_vnml.rgb = zeros(ny*ng,nx*ng,3,'uint8');
    imaged_data.seam_vnml.alp = zeros(ny*ng,nx*ng,1,'uint8');
    imaged_data.seam_vclr.rgb = zeros(ny*ng,nx*ng,3,'uint8');
    imaged_data.seam_vclr.alp = zeros(ny*ng,nx*ng,1,'uint8');
    imaged_data.seam_vtex.rgb = zeros(ny*32+1,nx*32+1,3,'uint8');
    imaged_data.seam_vtex.alp = zeros(ny*32+1,nx*32+1,1,'uint8');
    imaged_data.cmap.rgb = zeros(ny*256,nx*256,3,'uint8');
    imaged_data.cmap.alp = zeros(ny*256,nx*256,1,'uint8');
    
    % Ensure 'raw' export directory exists
    if(exist([opts.xprt_dir '\raw'],'dir')~=7)
        mkdir([opts.xprt_dir '\raw']);
    end
    
    % Scaling factor to convert heightmap to yards
    % (This many heightmap increments equals 1 yard)
    zfac = (cnst.cell_ng/cnst.cell_ny)/cnst.hinc_ng;
    
    % Scaling factor to convert heightmap to equal-to-horizontal units
    hfac = (cnst.cell_ny / cnst.cell_nh);
    
    % Scale height
    fhgt = ceil(mapped_data.fhgt./zfac);
    
    % Slope in units of rise-over-run, then get magnitude and direction
    [dzdx,dzdy] = gradient( mapped_data.fhgt./zfac );
    dz = dzdx + sqrt(-1).*dzdy;
    dz = dz./hfac;
    dz_mag = abs(dz); 
    dz_dir = angle(dz);
   
    % Scaling factor (trial-and-error derived...)
    % NOTE: Make pi/4 term settable?
    zshd = (dz_mag.*cos(dz_dir + pi/4));
    
    % Scale scaling factor
    %zshd = 128 + 128.*(zshd./max(abs(zshd(:))));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.01] BASEMAP
    imfile = [opts.xprt_dir '\raw\basemap.png'];
    if(exist(imfile,'file')~=2)
    
        % Easy since there's only 2 colors
        bmap_ind = uint8( mapped_data.fhgt(end:-1:1,:)>0 );
       	bmap_map = double([opts.gmap_wclr;opts.gmap_lclr])./255;
        imaged_data.bmap.rgb = ind2rgb(bmap_ind,bmap_map);
        imaged_data.bmap.rgb = uint8( imaged_data.bmap.rgb.*255 );
        imaged_data.bmap.alp = 255.*ones(ny*64,nx*64,'uint8');
        imwrite(imaged_data.bmap.rgb,imfile,'png','Alpha',imaged_data.bmap.alp);
    
    else
        [imaged_data.bmap.rgb,~,imaged_data.bmap.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.02] MINIMAP
    imfile = [opts.xprt_dir '\raw\minimap.png'];
    if(exist(imfile,'file')~=2)
    
        % This one is nice and easy
        imaged_data.mmap.rgb = uint8(repmat(mapped_data.wnam(end:-1:1,:),[1 1 3]));
        imaged_data.mmap.alp = 255.*ones(ny*9,nx*9,'uint8');
        imwrite(imaged_data.mmap.rgb,imfile,'png','Alpha',imaged_data.mmap.alp);
        
    else
        [imaged_data.mmap.rgb,~,imaged_data.mmap.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.03] REGION COLORS
    imfile = [opts.xprt_dir '\raw\region_colors.png'];
    if(isfield(merged_data.regn,'cnam'))
    if(exist(imfile,'file')~=2)
    
        % This one is easy...ish
        regn_ind = mapped_data.regn(end:-1:1,:)+1;
        regn_map = [0 0 0]';
        for i=1:numel(merged_data.regn)
            regn_map(:,end+1) = merged_data.regn(i).cnam(:,end);
        end
        imaged_data.regn.rgb = ind2rgb(regn_ind,single(regn_map')./255);
        imaged_data.regn.rgb = uint8( imaged_data.regn.rgb .* 255 );
        imaged_data.regn.alp = 255.*uint8( regn_ind~=1 );
        imwrite(imaged_data.regn.rgb,imfile,'png','Alpha',imaged_data.regn.alp);
        
    else
        [imaged_data.regn.rgb, ~, imaged_data.regn.alp] = imread(imfile);
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.04] VTEX INDICES
    imfile = [opts.xprt_dir '\raw\texture_indices.png'];
    if(exist(imfile,'file')~=2)
        
        % Number of indices
        n = max(mapped_data.vtex(:))+1;

        % Index colors
        clrs = tes3matlab_colormap_lines(double(n));
        clrs = clrs.*255;
        
        % Color indices
        imaged_data.vtex.rgb = ind2rgb(mapped_data.vtex(end:-1:1,:)+1,clrs);
        imaged_data.vtex.alp = 255.*uint8( mapped_data.vtex(end:-1:1,:)~=0 );
        
        % Write
        imwrite(imaged_data.vtex.rgb,imfile,'png','Alpha',imaged_data.vtex.alp);

    else
        [imaged_data.vtex.rgb,~,imaged_data.vtex.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.05] LTEX COLORS
    imfile = [opts.xprt_dir '\raw\texture_colors.png'];
    if(exist(imfile,'file')~=2)
    
        % Get list of colors from ltex file
        % NOTE: ASSUMES FIRST ENTRY IS _land_default
        fid = fopen(opts.lclr_file,'rt');
        tmp = textscan(fid,'%s %d %d %d');
        ltex_data = tmp{1};
        ltex_clr  = [tmp{2}, tmp{3}, tmp{4}];
        fclose(fid);
        
        % Link to intv's
        ii = zeros(numel(merged_data.ltex),1);
        for i=1:numel(merged_data.ltex)
            tmp = find( strcmp(ltex_data{i},{merged_data.ltex.data})==1 );
            if(~isempty(tmp))
                ii(i) = tmp;
            end
        end
        ii = ii+1;
        
        % Draw image
        ltex_ind = mapped_data.vtex(end:-1:1,:);
        ltex_map = single(ltex_clr(ii,:))./255;
        imaged_data.ltex.rgb = ind2rgb(ltex_ind,ltex_map);
        imaged_data.ltex.rgb = uint8( imaged_data.ltex.rgb.*255 );
        imaged_data.ltex.alp = 255.*ones(ny*16,nx*16,'uint8');
        
        % Save
        imwrite(imaged_data.ltex.rgb,imfile,'png','Alpha',imaged_data.ltex.alp);
        
    else
        [imaged_data.ltex.rgb,~,imaged_data.ltex.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.06] VERTEX COLORS
    imfile = [opts.xprt_dir '\raw\vertex_colors.png'];
    if(exist(imfile,'file')~=2)
    
        % Remove 65th rows & cols
        vclr = mapped_data.vclr;
        vclr(65:65:end,:,:) = [];
        vclr(:,65:65:end,:) = [];
        imaged_data.vclr.rgb = vclr(end:-1:1,:,:);
        
        % Set missing points to white
        tmp = sum(double(imaged_data.vclr.rgb),3);
        [ii,jj] = find(tmp==0);
        for i=1:numel(ii)
            imaged_data.vclr.rgb(ii(i),jj(i),:) = 255;
        end
        
        % Include alpha?
        tmp = sum(double(imaged_data.vclr.rgb),3);
        tmp(tmp==0 | tmp>=255*3) = 0;
        tmp(tmp~=0)=1;
        imaged_data.vclr.alp = tmp;
        
        % Write image
        imwrite(imaged_data.vclr.rgb,imfile,'png','Alpha',imaged_data.vclr.alp);
        
    else
        [imaged_data.vclr.rgb,~,imaged_data.vclr.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.07] VERTEX NORMALS
    imfile = [opts.xprt_dir '\raw\vertex_normals.png'];
    if(exist(imfile,'file')~=2)
    
        % This one's a little odd...
        vnml = mapped_data.vnml;
        vnml(65:65:end,:,:) = [];
        vnml(:,65:65:end,:) = [];
        imaged_data.vnml.rgb = uint8(single(vnml(end:-1:1,:,:))+127);
        imaged_data.vnml.alp = 255.*ones(ny*64,nx*64,'uint8');
        imwrite(imaged_data.vnml.rgb,imfile,'png','Alpha',imaged_data.vnml.alp);
        
    else
        [imaged_data.vnml.rgb,~,imaged_data.vnml.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [1.08] FACE HEIGHTS
    imfile = [opts.xprt_dir '\raw\face_heights.png'];
    if(exist(imfile,'file')~=2)
    
        % This one uses a fancy colormap
        fhgt_ind = fhgt(end:-1:1,:)+296;
        tmp      = tes3matlab_colormap_hgt;
        fhgt_map = interp1(tmp(:,1),tmp(:,2:4),tmp(end,1):1:tmp(1,1),'linear');
        imaged_data.fhgt.rgb = ind2rgb(fhgt_ind,fhgt_map./255);
        imaged_data.fhgt.alp = 255.*ones(ny*64,nx*64,'uint8');
        imwrite(imaged_data.fhgt.rgb,imfile,'png','Alpha',imaged_data.fhgt.alp);
        
    else
        [imaged_data.fhgt.rgb,~,imaged_data.fhgt.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.01] TOPOGRAPHY (SHADER)
    imfile = [opts.xprt_dir '\raw\topo_shade.png'];
    if(exist(imfile,'file')~=2)
        
        % Basemap with slopes shaded
        tshd = 128 + 128.*(zshd(end:-1:1,:)./10); %max(abs(zshd(:))));
        imaged_data.tshd.rgb = uint8(repmat(tshd,[1 1 3]));
        imaged_data.tshd.alp = 255.*ones(size(tshd),'uint8');
        imwrite(imaged_data.tshd.rgb,imfile,'png','Alpha',imaged_data.tshd.alp);
    
    else
        [imaged_data.tshd.rgb,~,imaged_data.tshd.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.02] TOPOGRAPHY (CONTOURS)
    imfile = [opts.xprt_dir '\raw\topo_contours.png'];
    if(exist(imfile,'file')~=2)
    
        % 10m and 100m heights
        z010m = -256.*ones(ny*64+2,nx*64+2);	z010m(2:end-1,2:end-1) = ceil(fhgt./010);
        z100m = -256.*ones(ny*64+2,nx*64+2);	z100m(2:end-1,2:end-1) = ceil(fhgt./100);
        
        % Contour image
        imaged_data.tcnt.rgb = zeros(ny*64,nx*64,'uint8');
        for iy=1:ny*64
        for ix=1:nx*64
            
            % This point and the surrounding points
            zz010 = z010m(iy+[0 1 2],ix+[0 1 2]);
            zz100 = z100m(iy+[0 1 2],ix+[0 1 2]);
            
            % Draw 10m contour?
            if(any(zz010(:)<z010m(iy+1,ix+1)))
                imaged_data.tcnt.rgb(iy,ix) = 64;
            end
            
            % Draw 100m contour?
            if(any(zz100(:)~=z100m(iy+1,ix+1)))
                imaged_data.tcnt.rgb(iy,ix) = 32;
            end
              
        end
        end
        
        % Fix
        imaged_data.tcnt.rgb = imaged_data.tcnt.rgb(end:-1:1,:);
        
        % Image alpha?
        imaged_data.tcnt.alp = 255.*uint8( imaged_data.tcnt.rgb~=0 );
        
        % Draw grayscale image directly?
        imwrite(imaged_data.tcnt.rgb,imfile,'png','Alpha',imaged_data.tcnt.alp);
        
    else
        [imaged_data.tcnt.rgb,~,imaged_data.tcnt.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.03] STEEP SLOPES
    imfile = [opts.xprt_dir '\raw\steep_slopes.png'];
    if(exist(imfile,'file')~=2)
        
        % Shade where slopes are >45deg (too steep to walk up)
        step = single(dz_mag(end:-1:1,:)<=1);
        imaged_data.step.rgb = uint8(255.*repmat(step,[1 1 3]));
        imaged_data.step.alp = 255.*uint8(step==0);
        imwrite(imaged_data.step.rgb,imfile,'png','Alpha',imaged_data.step.alp);
        
    else
        [imaged_data.step.rgb,~,imaged_data.step.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.04] ROADS & ABOVE-0 WATERS
    imfile = [opts.xprt_dir '\raw\roads.png'];
    if(exist(imfile,'file')~=2)
    
        % Get list of colors for roads
        fid = fopen(opts.rclr_file,'rt');
        tmp = textscan(fid,'%s %d %d %d');
        fclose(fid);
        road_data = tmp{1};
        road_clr  = [tmp{2}, tmp{3}, tmp{4}];
        
        % Link to intv's
        %ii = 0:numel(merged_data.ltex);
        jj = zeros(1,1);
        cc = zeros(1,3,'uint8');
        %ii = zeros(numel(road_data),2);
        for i=1:numel(road_data)
            tmp = find( strcmp(road_data{i},{merged_data.ltex.data})==1 );
            if(~isempty(tmp))
                jj(end+1)   = tmp;
                cc(end+1,:) = road_clr(i,:);
            end
        end
        tmp = mapped_data.vtex(end:-1:1,:);
        tmp(~ismember(tmp,jj)) = 1;
        tmp(mapped_data.vtex(end:-1:1,:)==0)=1;
        road_ind = zeros(ny*16,nx*16,'uint16');
        for i=1:numel(jj)
            [aa,bb] = find(tmp==jj(i));
            for j=1:numel(aa)
                road_ind(aa(j),bb(j)) = i;
            end
        end
        
        % Draw image
        imaged_data.road.rgb = ind2rgb(road_ind,double(cc)./255);
        imaged_data.road.rgb = uint8( imaged_data.road.rgb .* 255 );
        imaged_data.road.alp = 255.*uint8( road_ind~=0 );
        
        % Save
        imwrite(imaged_data.road.rgb,imfile,'png','Alpha',imaged_data.road.alp);
       
    else
        [imaged_data.road.rgb,~,imaged_data.road.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.05] VERTEX HEIGHT SEAMS
    imfile = [opts.xprt_dir '\raw\seams_vhgt.png'];
    if(exist(imfile,'file')~=2)
        
        % Seam image
        ng       = opts.seam_grid;
        seam_map = zeros(ny*ng,nx*ng,'uint8');
        
        % Loop through cells
        for iy=1:ny
        for ix=1:nx
            
            % Indices in vhgt
            ixv = (ix-1)*65+1 : ix*65;
            iyv = (iy-1)*65+1 : iy*65;
            
            % Indices in map
            ixm = (ix-1)*ng+1 : ix*ng;
            iym = (iy-1)*ng+1 : iy*ng;
            
            % Only do for cells with valid data
            tmp = mapped_data.vhgt(iyv,ixv);
            if(~all(tmp(:)==-256))
            
                % Left border
                if(ix==1)
                    brdr_lft = [-256.*ones(65,1), mapped_data.vhgt(iyv,ixv(1))];
                else
                    brdr_lft = mapped_data.vhgt(iyv,ixv(1)+[-1 0]);
                end
                diff_lft = abs(brdr_lft(:,2)-brdr_lft(:,1));
                if(~all(diff_lft==0))
                    seam_map(iym,ixm(1)) = uint8( max(diff_lft(:)) );
                end
                
                % Right border
                if(ix==nx)
                    brdr_rgt = [mapped_data.vhgt(iyv,ixv(end)), -256.*ones(65,1)];
                else
                    brdr_rgt = mapped_data.vhgt(iyv,ixv(end)+[0 1]);
                end
                diff_rgt = abs(brdr_rgt(:,2)-brdr_rgt(:,1));
                if(~all(diff_rgt==0))
                    seam_map(iym,ixm(end)) = uint8(max(diff_rgt(:)));
                end
                
                % Top border
                if(iy==1)
                    brdr_top = [-256.*ones(65,1), mapped_data.vhgt(iyv(1),ixv)'];
                else
                    brdr_top = mapped_data.vhgt(iyv(1)+[-1 0],ixv)';
                end
                diff_top = abs(brdr_top(:,2)-brdr_top(:,1));
                if(~all(diff_top==0))
                    seam_map(iym(1),ixm) = uint8(max(diff_top(:)));
                end

                % Bottom border
                if(iy==ny)
                    brdr_btm = [mapped_data.vhgt(iyv(end),ixv)', -256.*ones(65,1)];
                else
                    brdr_btm = mapped_data.vhgt(iyv(end)+[0 1],ixv)';
                end
                diff_btm = abs(brdr_btm(:,2)-brdr_btm(:,1));
                if(~all(diff_btm==0))
                    seam_map(iym(end),ixm) = uint8(max(diff_btm(:)));
                end
                
            end
            
        end
        end
        
        % Fix
        seam_map = seam_map(end:-1:1,:);
        
        % Write?
        imaged_data.seam_vhgt.rgb(:,:,1) = seam_map;
        %imaged_data.seam_vhgt.rgb(:,:,2) = 0;
        %imaged_data.seam_vhgt.rgb(:,:,3) = 0;
        imaged_data.seam_vhgt.alp = 255.*uint8(seam_map~=0);
        imwrite(imaged_data.seam_vhgt.rgb,imfile,'Alpha',imaged_data.seam_vhgt.alp);
        
    else
        [imaged_data.seam_vhgt.rgb,~,imaged_data.seam_vhgt.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.06] VERTEX COLOR SEAMS
    imfile = [opts.xprt_dir '\raw\seams_vclr.png'];
    if(exist(imfile,'file')~=2)
        
        % Seam image
        ng       = opts.seam_grid;
        seam_map = zeros(ny*ng,nx*ng,'uint8');
        
        % vclr
        tmp = sum(single(mapped_data.vclr(end:-1:1,:,:)),3);
        
        % Loop through cells
        for iy=1:ny
        for ix=1:nx
            
            % Indices in vhgt
            ixv = (ix-1)*65+1 : ix*65;
            iyv = (iy-1)*65+1 : iy*65;
            
            % Indices in map
            ixm = (ix-1)*ng+1 : ix*ng;
            iym = (iy-1)*ng+1 : iy*ng;
            
            % Only do for cells with valid data
            tmpp = tmp(iyv,ixv);
            if(~all(tmpp(:)==0))
            
                % Left border
                if(ix==1)
                    brdr_lft = [765.*ones(65,1), tmp(iyv,ixv(1))];
                else
                    brdr_lft = tmp(iyv,ixv(1)+[-1 0]);
                end
                diff_lft = abs(brdr_lft(:,2)-brdr_lft(:,1));
                if(~all(diff_lft==0))
                    seam_map(iym,ixm(1)) = uint8( max(diff_lft(:)) );
                end
                
                % Right border
                if(ix==nx)
                    brdr_rgt = [tmp(iyv,ixv(end)), 765.*ones(65,1)];
                else
                    brdr_rgt = tmp(iyv,ixv(end)+[0 1]);
                end
                diff_rgt = abs(brdr_rgt(:,2)-brdr_rgt(:,1));
                if(~all(diff_rgt==0))
                    seam_map(iym,ixm(end)) = uint8(max(diff_rgt(:)));
                end
                
                % Top border
                if(iy==1)
                    brdr_top = [765.*ones(65,1), tmp(iyv(1),ixv)'];
                else
                    brdr_top = tmp(iyv(1)+[-1 0],ixv)';
                end
                diff_top = abs(brdr_top(:,2)-brdr_top(:,1));
                if(~all(diff_top==0))
                    seam_map(iym(1),ixm) = uint8(max(diff_top(:)));
                end

                % Bottom border
                if(iy==ny)
                    brdr_btm = [tmp(iyv(end),ixv)', 765.*ones(65,1)];
                else
                    brdr_btm = tmp(iyv(end)+[0 1],ixv)';
                end
                diff_btm = abs(brdr_btm(:,2)-brdr_btm(:,1));
                if(~all(diff_btm==0))
                    seam_map(iym(end),ixm) = uint8(max(diff_btm(:)));
                end
                
            end
            
        end
        end
        
        % Write?
        imaged_data.seam_vclr.rgb(:,:,1) = seam_map;
        %imaged_data.seam_vclr(:,:,end+1) = 0;
        %imaged_data.seam_vclr(:,:,end+1) = 0;
        imaged_data.seam_vclr.alp = 255.*uint8(seam_map~=0);
        imwrite(imaged_data.seam_vclr.rgb,imfile,'Alpha',imaged_data.seam_vclr.alp);
        
    else
        [imaged_data.seam_vclr.rgb,~,imaged_data.seam_vclr.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.07] LAND TEXTURE VERTEX SEAMS (>3 at a vertex)
    imfile = [opts.xprt_dir '\raw\seams_vtex.png'];
    if(exist(imfile,'file')~=2)
    
        % Seam image
        ng       = opts.seam_grid;
        seam_map = zeros(ny*32+1,nx*32+1,'uint8');
        
        % VTEX
        tmp = mapped_data.vtex(end:-1:1,:);
        
        % Loop
        for iy=1:size(tmp,1)-1
        for ix=1:size(tmp,2)-1
        
            % 4 surrounding vertices
            tmpp = tmp(iy:(iy+1),ix:(ix+1));
            n    = numel(unique(tmpp));
            
            % If 3 or 4, save to image
            if(n>2)
                iiy = iy*2;
                iix = ix*2;
                seam_map(iiy,iix) = n-2;
            end
            
        end
        end
        
        % Draw image
        cc = uint8([0 0 0; 64 0 0; 255 0 0]);
        imaged_data.seam_vtex.rgb = ind2rgb(seam_map,cc)./255;
        imaged_data.seam_vtex.alp = 255.*uint8(seam_map~=0);
        
        % Save
        imwrite(imaged_data.seam_vtex.rgb,imfile,'png','Alpha',imaged_data.seam_vtex.alp);
            
    else
        [imaged_data.seam_vtex.rgb,~,imaged_data.seam_vtex.alp] = imread(imfile);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % [2.08] CREATEMAPS MERGED
%     imfile = [opts.xprt_dir '\raw\createmaps_merged.png'];
%     if(exist(imfile,'file')~=2)
%     if(~isempty(opts.cmap_dir))
%      
%         % Possible x's and y's
%         x = [merged_data.land.x];   x = min(x):1:max(x);
%         y = [merged_data.land.y];   y = min(y):1:max(y);
%         
%         % Get CreateMaps images within this x and y range
%         cmap_data = tes3matlab_get_cmap([merged_data.land.x;merged_data.land.y],opts.cmap_dir);
%         xy = [cmap_data.x,;cmap_data.y];
%         
%         % Grid
%         ng = opts.cmap_size;
% 
%         % Loop
%         for iy=1:ny
%         for ix=1:nx
%             
%             % Continue if image exists for this cell
%             ii = find(xy(1,:)==x(ix) & xy(2,:)==y(iy));
%             if(~isempty(ii))
%             if(~isempty(cmap_data(ii).file))
%             
%                 % Add
%                 iiy = (iy-1)*ng+1 : iy*ng;
%                 iix = (ix-1)*ng+1 : ix*ng;
%                 img = imread([opts.cmap_dir '\' cmap_data(ii).file]);
%                 imaged_data.cmap.rgb(iiy,iix,:) = img(end:-1:1,:,:);
%                 imaged_data.cmap.alp(iiy,iix) = 255;
%                 
%             end
%             end
%         
%         end
%         end
%         imaged_data.cmap.rgb = imaged_data.cmap.rgb(end:-1:1,:,:);
%         imaged_data.cmap.alp = imaged_data.cmap.alp(end:-1:1,:);
%         
%         % Save?
%         imwrite(imaged_data.cmap.rgb,imfile,'png','Alpha',imaged_data.cmap.alp);
%         
%     end
%     else
%         [imaged_data.cmap.rgb,~,imaged_data.cmap.alp] = imread(imfile);
%     end
    
end