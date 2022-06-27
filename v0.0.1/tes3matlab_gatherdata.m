function file_data = tes3matlab_gatherdata(opts)
% file_data = tes3matlab_gather(opts)
% Gathers data from TES3 files listed in the provided text file and saves
% the data from each file in the project's data directory.  If the data was
% extracted previously, the data is loaded from the project data file.
% 
% If the original file gets updated and you wish to re-run the extraction, 
% make sure to first delete the associated data file in the project data 
% directory!
% 
% by Spider the spider

    %=====================================================================%

    % Try to open list_file
    fid = fopen(opts.list_file,'rt');
    if(fid==-1);	error(['Unable to open file ' list_file ' for reading.']);  end
    
    % Read in line by line
    file_list = textscan(fid,'%s\n','delimiter','');
    file_list = file_list{1};

    % Close file
    fclose(fid);
    clear fid;
    
    %=====================================================================%
    % Ensure files exist before attempting to read in their data
    n_file = numel(file_list);
    for i=1:n_file

        % Attempt to open the file.
        % If successful, ensure it's a TES3 file.
        fid = fopen(file_list{i},'rb');
        if(fid<=0)
            error(['Unable to open file ' file_list{i} ' for reading.']);
        else
            file_signature = native2unicode( fread(fid,4,'ubit8') )';
            if(strcmp(file_signature,'TES3')~=1)
                error(['Invalid file signature: ' file_signature ', this version only works for TES3 files!']);
            end
        end
        fclose(fid);

    end %for i=1:n_file
    clear i fid file_signature;
    %=====================================================================%

    %=====================================================================%
    % Loop through files and gather their data
    for i=1:n_file

        % Check if data has been read in before
        % If not, extract the raw data and save to a data file.
        file_name = strsplit(file_list{i},{'\','/','.'});
        file_name = file_name{end-1};
        data_file = [opts.dir_data '\' file_name '.mat'];   % NOTE: IN FUTURE WILL BE NETCDF FOR MORE UNIVERSAL ACCESS!
        if(exist(data_file,'file')~=2)
            tes3matlab_extractdata(file_list{i},data_file);
        else
            warning(['Data from ' file_list{i} ' has already been extracted previously and saved to file ' data_file '.']);
        end
 
    end %for i=1:n_file
    clear i file_name data_file;
    %=====================================================================%
    
end