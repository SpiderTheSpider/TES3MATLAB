function cmap = tes3matlab_get_cmap(xy,cmap_dir)

    % Number of locations
    x = xy(1,:);
    y = xy(2,:);
    n = numel(x);
    
    % cmap structure
    cmap = struct;
    
    % Loop
    for i=1:n

        % This position
        cmap(i).x = x(i);
        cmap(i).y = y(i);
        
        % Look for CreateMaps image with these coordinates
        test = ls([cmap_dir '\*(' num2str(x(i)) ',' num2str(y(i)) ')*']);
        if(~isempty(test))
            
            % File name and image
            cmap(i).file = test;
            cmap(i).img  = imread([cmap_dir '\' test]);
            
            % Location name
            ii = strfind(test,'(');
            cmap(i).name = lower(deblank( test(1:ii-1) ));
            
        end

    end

end