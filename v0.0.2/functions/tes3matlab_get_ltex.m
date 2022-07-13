function ltex = tes3matlab_get_ltex(ltex_list,ltex_dir,lclr_file)

    % Oops, add _land_default
    ltex_list{end+1} = '_land_default';

    % Output data structure
    ltex = struct;
    n    = numel(ltex_list);
    b    = 0;                   % For missing textures
    
    % Loop through entries
    for i=1:n
        
        % File name
        ltex(i).data = ltex_list{i};
        
        % Look for image
        test = ls([ltex_dir '\' ltex_list{i} '.*']);
        if(~isempty(test))
            ltex(i).img = imread([ltex_dir '\' test]);
            ltex(i).clr = immean(ltex(i).img);
        else
            ltex(i).img = [];
            ltex(i).clr = uint8([b b b]);
            b           = b + 2;
        end

    end
    
    % If list of colors doesn't yet exist, generate it
    if(exist(lclr_file,'file')~=2)
       fid = fopen(lclr_file,'wt');
       for i=1:numel(ltex)
           fprintf(fid,'%30s %0.3d %0.3d %0.3d\n',ltex(i).data,ltex(i).clr);
       end
       fclose(fid);
    end

end