function data = tes3matlab_loaddata(data_file,var_to_get)

    % Load data from data_file
    % NOTE: This simple code is a placeholder for a future version that
    % will save to a more accessible format (e.g., NetCDF)
    load_str = 'data = load(data_file';
    for i=1:numel(var_to_get)
        load_str = [load_str ',''' var_to_get{i} ''''];
    end
    load_str = [load_str ');'];
    eval(load_str);
    
end