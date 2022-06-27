function [cmap_img, cmap_info] = tes3matlab_get_cmap(xy,cmap_dir)

    % Locations
    n = size(xy,2);
    
    % CreateMaps data
    cmap_img = cell(n,1);
    cmap_info = cell(n,1);
    
    % Loop through xy's
    for i=1:n

        % See if CreateMaps images exist
        cmap_file = ls([cmap_dir '\*(' num2str(xy(1,i)) ',' num2str(xy(2,i)) ')*.bmp']);
        if(~isempty(cmap_file))
            cmap_info{i} = cmap_file;
            cmap_img{i}  = imread([cmap_dir '\' cmap_file]); 
        end
        
    end
    clear i n cmap_file;

end