%function tes3matlab_drawautomap(mapped_data,imaged_data,opts)
% 
% The Automap is a 64x64 pixels-per-grid image file of .tiff layers:
% (1) Basemap (land/sea mask)
% (2) Cell boundaries (shader to show cell boundaries)
% (3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % World dimensions
    [ny,nx] = size(imaged_data.regn.alp);
    ng      = 40;
    
    % Init the Automap tiff
    auto_file = [opts.proj_dir '\automap.tif'];
    if(exist(auto_file,'file')==2)
        delete(auto_file);
    end
    tf = Tiff(auto_file,'w');

    % Some initial tags
    ts = struct;
    ts.ImageLength          = ny*ng;
    ts.ImageWidth           = nx*ng;
    ts.Photometric          = Tiff.Photometric.RGB;
    ts.BitsPerSample        = 8;
    ts.SamplesPerPixel      = 4; % 3 + alpha channel
    ts.PlanarConfiguration	= Tiff.PlanarConfiguration.Chunky; 
	ts.Software             = 'MATLAB';
    %ts.TileLength           = ng;
    %ts.TileWidth            = ng;
    ts.ExtraSamples         = Tiff.ExtraSamples.AssociatedAlpha;
    %ts.Compression          = Tiff.Compression.None;
    %ts.JPEGColorMode        = Tiff.JPEGColorMode.RGB;
    %ts.JPEGQuality          = 100;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % First layer: the basemap
    a = imresize(imaged_data.bmap.rgb,[ny*ng nx*ng],'nearest');
    a(:,:,end+1) = 255.*ones(ny*ng,nx*ng,'uint8');
    %imwrite(a,auto_file);
    ts.PageName = 'Basemap';
    setTag(tf,ts);
    write(tf,a);
    
    % Second layer: cell defs
    a = 255.*ones(ny*ng,nx*ng,4,'uint8');
    for iy=1:ny
    for ix=1:nx
        iiy = (iy-1)*ng+1 : iy*ng;
        iix = (ix-1)*ng+1 : ix*ng;
        if(mod(iy+ix,2)==0)
            a(iiy,iix,1:3) = 0;
        else
            a(iiy,iix,1:3) = 255;
        end
    end
    end
    writeDirectory(tf); 
    ts.PageName = 'Cell_Shade';
    setTag(tf,ts);  
    write(tf,a);
    
    % ??? layer: heightmap shade
    a = imresize(imaged_data.tshd.rgb,[ny*ng nx*ng],'nearest');
    b = imresize(imaged_data.tshd.alp,[ny*ng nx*ng],'nearest');
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'Height_Overlay';
    setTag(tf,ts);  
    write(tf,a);
    
    % ??? layer: minimap
    a = imresize(imaged_data.mmap.rgb,[ny*ng nx*ng],'nearest');
    b = imresize(imaged_data.mmap.alp,[ny*ng nx*ng],'nearest');
    a(:,:,end+1) = b;
    %imwrite(a,auto_file,'WriteMode','append');
    writeDirectory(tf); 
    ts.PageName = 'Minimap';
    setTag(tf,ts);  
    write(tf,a);
    
%     % ??? layer: region colors
%     a = imresize(imaged_data.regn.rgb,[ny*ng nx*ng],'nearest');
%     b = imresize(imaged_data.regn.alp,[ny*ng nx*ng],'nearest');
%     a(:,:,end+1) = b;
%     %imwrite(a,auto_file,'WriteMode','append');
%     writeDirectory(tf); 
%     ts.PageName = 'Region_Colors';
%     setTag(tf,ts);  
%     write(tf,a);
    
%     % ??? layer: land texture colors
%     a = imresize(imaged_data.ltex.rgb,[ny*ng nx*ng],'nearest');
%     b = imresize(imaged_data.ltex.alp,[ny*ng nx*ng],'nearest');
%     a(:,:,end+1) = b;
%     writeDirectory(tf); 
%     ts.PageName = 'LandTexture_Colors';
%     setTag(tf,ts);  
%     write(tf,a);
    
    % ??? layer: heightmap colors
    a = imresize(imaged_data.fhgt.rgb,[ny*ng nx*ng],'nearest');
    b = imresize(imaged_data.fhgt.alp,[ny*ng nx*ng],'nearest');
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'Height_Colors';
    setTag(tf,ts);  
    write(tf,a);

    % ??? layer: heightmap contours
    a = imresize(imaged_data.tcnt.rgb,[ny*ng nx*ng],'nearest');
    b = imresize(imaged_data.tcnt.alp,[ny*ng nx*ng],'nearest');
    a(:,:,end+1) = a(:,:,1);
    a(:,:,end+1) = a(:,:,1);
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'Height_Contours';
    setTag(tf,ts);  
    write(tf,a);
    
%     % ??? layer: road highlights
%     a = imresize(imaged_data.road.rgb,[ny*ng nx*ng],'nearest');
%     b = imresize(imaged_data.road.alp,[ny*ng nx*ng],'nearest');
%     a(:,:,end+1) = b;
%     writeDirectory(tf); 
%     ts.PageName = 'Road_Highlights';
%     setTag(tf,ts);  
%     write(tf,a);
    
    % ??? layer: steeps
    a = imresize(imaged_data.step.rgb,[ny*ng nx*ng]);
    b = imresize(imaged_data.step.alp,[ny*ng nx*ng]);
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'Steepness_ID';
    setTag(tf,ts);  
    write(tf,a);
    
    % ??? layer: vertex colors
    a = imresize(imaged_data.vclr.rgb,[ny*ng nx*ng]);
    b = imresize(255.*imaged_data.vclr.alp,[ny*ng nx*ng]);
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'Vertex_Colors';
    setTag(tf,ts);  
    write(tf,a);

    % ??? layer: vertex normals
    a = imresize(imaged_data.vnml.rgb,[ny*ng nx*ng]);
    b = imresize(imaged_data.vnml.alp,[ny*ng nx*ng]);
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'Vertex_Normals';
    setTag(tf,ts);  
    write(tf,a);
    
%     % ??? layer: createmaps
%     a = imresize(imaged_data.cmap.rgb,[ny*ng nx*ng]);
%     b = imresize(imaged_data.cmap.alp,[ny*ng nx*ng]);
%     a(:,:,end+1) = b;
%     writeDirectory(tf); 
%     ts.PageName = 'CreateMaps_Merged';
%     setTag(tf,ts);  
%     write(tf,a);
    
    % ??? layer: seam test: heightmap
    a = imresize(imaged_data.seam_vhgt.rgb,[ny*ng nx*ng]);
    b = imresize(imaged_data.seam_vhgt.alp,[ny*ng nx*ng]);
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'SeamTest_VertexHeights';
    setTag(tf,ts);  
    write(tf,a);
    
    % ??? layer: seam test: vertex colors
    a = imresize(imaged_data.seam_vclr.rgb,[ny*ng nx*ng]);
    b = imresize(imaged_data.seam_vclr.alp,[ny*ng nx*ng]);
    a(:,:,end+1) = b;
    writeDirectory(tf); 
    ts.PageName = 'SeamTest_VertexColors';
    setTag(tf,ts);  
    write(tf,a);
    
    % ??? layer: seam test: vertex land textures
    % hmmm, how to resize this one...
    
    
    % Done
    close(tf);

%end