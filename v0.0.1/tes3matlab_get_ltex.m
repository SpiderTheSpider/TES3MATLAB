function [ltex_img, ltex_clr] = tes3matlab_get_ltex(ltex_list,ltex_dir)

    % Load
    ltex_list{end+1} = '_land_default';
    n_ltex           = numel(ltex_list);
    ltex_img         = cell(n_ltex,1);
    ltex_clr         = zeros(n_ltex,3,'uint8');
    tmp              = 4;
    for i=1:n_ltex
        this_file = [ltex_dir '\' ltex_list{i} '.png'];
        if(exist(this_file,'file')==2)
            raw_img = imread(this_file);
            ltex_img{i}   = raw_img;
            for j=1:3;	ltex_clr(i,j) = median(median(raw_img(:,:,j))); end
        else
            ltex_img{i}   = [];
            ltex_clr(i,:) = [1 1 1].*tmp;
            tmp = tmp + 2;
        end
        clear raw_img this_file;
    end
    clear i tmp;

end