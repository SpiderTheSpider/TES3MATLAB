function file_list = my_ls(file_dir)

    % NOTE: This works for windows-based machines, but not for UNIX-based
    % FUTURE: Add switch/case dependent on system architexture
    tmp = ls(file_dir);
    file_list = cell(size(tmp,1),1);
    for i=1:size(tmp,1)
        file_list{i} = deblank(tmp(i,:));
    end

end