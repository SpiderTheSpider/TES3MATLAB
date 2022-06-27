%function tes3matlab_draw_gridmap(merged_data,map_data)

    %=====================================================================%
    % INIT
    
    % Gridmap colors
    land_color = [222 199 176];
    watr_color = [133 164 179];

    % World size in cells
    xy = [merged_data.land.xy];
    x  = min(xy(1,:)) : max(xy(1,:));   nx = numel(x);  x0 = min(x);
    y  = min(xy(2,:)) : max(xy(2,:));   ny = numel(y);  y0 = min(y);
    
    % Cell size in pixels
    ng = opts.grid_size;
    
    % Meshgrid for interpolating other maps to gridmap's size
    [xg,yg] = meshgrid(0:(nx*ng-1),0:(ny*ng-1));    %xg=xg'; yg=yg';
    
    %=====================================================================%
    % PROCESSING
    
    % Interpolate heightmap to gridmap scale
    nh = opts.dims.cell_nh;
    if(nh~=ng)
        [xh,yh]   = meshgrid(linspace(0,nx*ng-1,nx*nh),linspace(0,ny*ng-1,ny*nh));
        grid_fhgt = interp2(xh,yh,map_data.fhgt,xg,yg,'linear');
    else
        grid_fhgt = map_data.fhgt;
    end
    
    % Heightmap in units of yards
    hscale = (opts.dims.cell_ng/opts.dims.hinc_ng)/opts.dims.cell_ny;
    grid_fhgt = grid_fhgt./hscale;
    
    % Land-sea mask
    grid_mask = grid_fhgt > 0;
    
    % Height contours to 1yd, 10yd
    grid_cntr_01yd = ceil(grid_fhgt);
    grid_cntr_10yd = ceil(grid_fhgt./10);
    
    % Compute gradients
    [dzdx,dzdy] = gradient(grid_fhgt);  % Slopes in terms of derivatives dz/dx, dz/dy
    dz = dzdx + sqrt(-1).*dzdy;         % Vector form
    
    % Convert gradients to units of rise/run
    dz = dz*(ng/opts.dims.cell_ny);
    
    % Magnitude and direction as angle counter-clockwise from Eastward
    dz_mag = abs(dz);                   % Magnitude
    dz_dir = angle(dz);                 % Direction (NOTE: Visually, -dz_dir (cw from Westward) looks better!)
    
    % Land texture map
    nt        = opts.dims.cell_nt;
    if(nt~=ng)
        [xt,yt]   = meshgrid(linspace(0,nx*ng-1,nx*nt),linspace(0,ny*ng-1,ny*nt));
        grid_tclr = interp2(xt,yt,map_data.ltex,xg,yg,'nearest');
    else
        grid_tclr = map_data.ltex;
    end
    
    %=====================================================================%
    % IMAGING
    
    % Land texture colors
    img_tclr = zeros(ny*ng,nx*ng,3,'uint8');
    for i=1:ny*ng
    for j=1:nx*ng
        
        % This LTEX index
        ii = grid_tclr(i,j);
        if(ii==0); ii = size(ltex_clr,1);   end
        
        % Draw mean land texture color
        img_tclr(i,j,:) = ltex_clr(ii,:);
        
    end
    end
    clear i j ii;
    
    % Player cannot ascend slopes greater than ~46 (use 45 as estimate
    % since it's easy to find; where rise = run)
    grid_slop = dz_mag<1;
    img_slop  = uint8(repmat(grid_slop,[1 1 3])).*255;
    img_slos  = uint8(repmat(atand(dz)./90,[1 1 3]).*255);
    %imshow(grid_slop);
    %imshow(img_tclr.*repmat(uint8(grid_slop),[1 1 3])); set(gca,'ydir','normal');
    
    % Vertex colors
    tmp = map_data.vclr;    tmp(65:65:end,:,:) = [];    tmp(:,65:65:end,:) = [];
    if(nh~=ng)
        img_vclr = zeros(ny*ng,nx*ng,3,'uint8');
        for i=1:3
            img_vclr(:,:,i) = uint8(interp2(xh,yh,single(tmp(:,:,i)),xg,yg,'linear'));
        end
    else
        img_vclr = tmp;
    end
    clear tmp i;
    
    % Landscape shading (note: computed as vclr * landscape)
    tmp = (double(img_tclr)./255) .* (double(img_vclr)./255);
    img_shad = uint8(tmp.*255);
    
    % Landscape mask
    img_mask = zeros(ny*ng,nx*ng,3,'uint8');
    for i=1:ny*ng
    for j=1:nx*ng
    if(grid_fhgt(i,j,:)>0)
        img_mask(i,j,:) = land_color;
    else
        img_mask(i,j,:) = watr_color;
    end
    end
    end
    
    % Region colors
    img_rclr = zeros(ny*ng,nx*ng,3,'uint8');
    for i=1:ny*ng
    for j=1:nx*ng
        img_rclr(i,j,:) = map_data.rclr(ceil(i/ng),ceil(j/ng),:);
    end
    end
    
    % CreateMaps
    nc = opts.dims.cell_nc;
    img_cmap = zeros(ny*ng,nx*ng,3,'uint8');
    for iy=1:ny
    for ix=1:nx
    	ii = find(xy(1,:)==x(ix) & xy(2,:)==y(iy));
        if(~isempty(ii))
        if(~isempty(cmap_img{ii}))
            xx = (ix-1)*ng+1 : ix*ng;
            yy = (iy-1)*ng+1 : iy*ng;
            if(nc~=ng)
                img_cmap(yy,xx,:) = imresize(cmap_img{ii}(end:-1:1,:,:),[ng ng]);
            else
                img_cmap(yy,xx,:) = cmap_img{ii}(end:-1:1,:,:);
            end
        end
        end
    end
    end
    clear iy ix ii xx yy;

    % Landscape height by color, line contours, and bevel contours
    z        = floor(grid_fhgt); z(z<-200)=-200; z(z>299)=299;
    z0       = -201;
    zclrmap  = demcmap([-200 300],500);
    img_zclr = zeros(ny*ng,nx*ng,3,'uint8');
    img_c01y = ones(ny*ng,nx*ng,3,'uint8').*255;
    img_c10y = ones(ny*ng,nx*ng,3,'uint8').*255;
    for iy=1:ny*ng
    for ix=1:nx*ng
        
    	% Height as color
        img_zclr(iy,ix,:) = uint8(zclrmap(z(iy,ix)-z0,:).*255);
        
%         % Height as 1yd contour
%         this_z = grid_cntr_01yd(iy,ix);
%         next_z = []; 
%         if(iy>1);     next_z = [next_z, grid_cntr_01yd(iy-1,ix)]; end
%         if(iy<ny*ng); next_z = [next_z, grid_cntr_01yd(iy+1,ix)]; end
%         if(ix>1);     next_z = [next_z, grid_cntr_01yd(iy,ix-1)]; end
%         if(ix<nx*ng); next_z = [next_z, grid_cntr_01yd(iy,ix+1)]; end
%         if(any(next_z<this_z))
%             img_c01y(iy,ix,:) = 128;
%         end
        
        % Height as 10yd contour
        this_z = grid_cntr_10yd(iy,ix);
        next_z = []; 
        if(iy>1);     next_z = [next_z, grid_cntr_10yd(iy-1,ix)]; end
        if(iy<ny*ng); next_z = [next_z, grid_cntr_10yd(iy+1,ix)]; end
        if(ix>1);     next_z = [next_z, grid_cntr_10yd(iy,ix-1)]; end
        if(ix<nx*ng); next_z = [next_z, grid_cntr_10yd(iy,ix+1)]; end
        if(any(next_z<this_z))
            img_c10y(iy,ix,:) = 64;
        end
        
    end
    end
    
    % Sloped shading using dz_dir
    % shading = mag*sin(dir)
    tmp = 1+(dz_mag.*cos(dz_dir+pi/4))./4;
    tmp = (double(img_mask)./255) .* tmp;
    img_slor = uint8(tmp.*255);
    %imshow(img_slor(end:-1:1,:,:))
    
    % Region colors
    fig = figure('visible','off');
    axx = axes; hold on;
    rr  = unique(map_data.regn(:));
    rr  = rr(rr>0);
    for i=1:numel(merged_data.regn)
        scatter(NaN,NaN,24,'s','filled');
    end; hold off;
    [~,li] = legend({merged_data.regn(rr).fnam});
    li = findobj(li,'Type','line');
    li = findobj(li,'Marker','none','-xor');
    set(li,'MarkerSize',20);
    set(fig,'visible','on');


    
    % Test look
    %cmapsea  = [0 0 0; 0 1/8 1/4; 0 1/4 1/2; 0 3/8 3/4; 0 1/2 1];
    %cmapland = [24 48 28; 208 224 255; 64 24 24]./255;
    %clrmap   = demcmap([-200 350],450,cmapsea,cmapland);
    %zz = ceil(grid_fhgt(end:-1:1,:)); zr = [max(zz(:)) min(zz(:))];
    %imshow(zz,[]); colormap(demcmap(zr)); colorbar;
    
%     % Region shading (note: computed as mask*region color)
%     tmp = (double(img_mask)./255) .* (double(img_rclr)./255);
%     img_shar = uint8(tmp.*255);
%     imshow(img_shar(end:-1:1,:,:));

%     % Shaded land textures times land-sea mask
%     tmp = (double(img_shad)./255) .* (double(img_mask)./255);
%     img_test = uint8(tmp.*255);
%     imshow(img_test(end:-1:1,:,:));

    % Export
    xprt_dir = [opts.dir_xprt '\gridmap'];
    if(exist(xprt_dir,'dir')~=7);   mkdir(xprt_dir);    end
    imwrite(img_mask(end:-1:1,:,:),[xprt_dir '\landmass_flat.png'],'png');
    imwrite(img_rclr(end:-1:1,:,:),[xprt_dir '\regions.png'],'png');
    imwrite(img_tclr(end:-1:1,:,:),[xprt_dir '\ltex_solid.png'],'png');
    imwrite(img_shad(end:-1:1,:,:),[xprt_dir '\ltex_shade.png'],'png');
    imwrite(img_vclr(end:-1:1,:,:),[xprt_dir '\vertex_colors.png'],'png');
    imwrite(img_slop(end:-1:1,:,:),[xprt_dir '\slopes_solid.png'],'png');
    imwrite(img_slor(end:-1:1,:,:),[xprt_dir '\landmass_shaped.png'],'png');
    imwrite(img_slos(end:-1:1,:,:),[xprt_dir '\slopes_shade.png'],'png');
    imwrite(img_cmap(end:-1:1,:,:),[xprt_dir '\createmaps.png'],'png');
    imwrite(img_zclr(end:-1:1,:,:),[xprt_dir '\height_color.png'],'png');
    imwrite(img_c01y(end:-1:1,:,:),[xprt_dir '\height_contour_01yd.png'],'png');
    imwrite(img_c10y(end:-1:1,:,:),[xprt_dir '\height_contour_10yd.png'],'png');
    