function tes3matlab_savedata(data,save_file) %,opts)

    % Save data to save_file
    % NOTE: This simple code is a placeholder for a future version that
    % will save to a more accessible format (e.g., NetCDF)
    var_to_save = fieldnames(data);
    save_str    = 'save(save_file';
    for i=1:numel(var_to_save)
        eval([var_to_save{i} '=data.' var_to_save{i} ';']);
        save_str = [save_str ',''' var_to_save{i} ''''];
    end
    save_str = [save_str ');'];
    eval(save_str);

end