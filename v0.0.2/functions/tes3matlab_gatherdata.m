function file_data = tes3matlab_gatherdata(opts,cnst)
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
    % INIT
    
    % Init variables
    fid               = [];         % File identifier
    listof_tes3files = {};          % Cell list of paths to tes3 files
    listof_datafiles  = {};         % Cell list of paths to data files
    n_files           = uint8(0);   % Number of input files
    file_signature    = '';         % File signature, should be TES3
    tes3_file         = '';         % Current tes3 file
    data_file         = '';         % Current data file
    file_data         = struct;
     file_data.cell   = struct;
     file_data.land   = struct;
     file_data.ltex   = struct;
     file_data.regn   = struct;
    
    %=====================================================================%

    % Try to open list_file
    fid = fopen(opts.proj_infile,'rt');
    if(fid==-1);	error(['Unable to open file ' opts.proj_infile ' for reading.']);  end
    
    % Read in line by line
    listof_tes3files = textscan(fid,'%s\n','delimiter','');
    if(numel(listof_tes3files)==1);	listof_tes3files = listof_tes3files{1};	end
    if(~iscell(listof_tes3files));  listof_tes3files = {listof_tes3files};  end

    % Close file
    fclose(fid);
    
    %=====================================================================%
    % Ensure files exist before attempting to read in their data
    n_files(1)       = numel(listof_tes3files);
    listof_datafiles = cell(n_files,1);
    for i=1:n_files

        % Attempt to open the file.
        % If successful, ensure it's a TES3 file.
        fid = fopen(listof_tes3files{i},'rb');
        if(fid<=0)
            error(['Unable to open file ' listof_tes3files{i} ' for reading.']);
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
    % Generate list of data files associated with this project
    listof_datafiles = cell(n_files,1);
    for i=1:n_files

        % Get this file's name and use it to generate the data file name
        tes3_file = strsplit(listof_tes3files{i},{'\','/','.'});
        tes3_file = tes3_file{end-1};
        data_file = [opts.data_dir '\' tes3_file '.mat'];   % NOTE: IN FUTURE WILL BE NETCDF FOR MORE UNIVERSAL ACCESS!
        
        % Save data file name
        listof_datafiles{i} = data_file;
        
        % If the data file doesn't exist, extract data from the input file
        % Else, load the data
        if(exist(data_file,'file')~=2)
            this_data = tes3matlab_extractdata(listof_tes3files{i},opts,cnst);
        else
            this_data = tes3matlab_loaddata(data_file,{'cell','land','ltex','regn','file'});
        end
        
        % Save to the file_data structure
        file_data(i).file = this_data.file;
        file_data(i).cell = this_data.cell;
        file_data(i).land = this_data.land;
        file_data(i).ltex = this_data.ltex;
        file_data(i).regn = this_data.regn;    

    end %for i=1:n_files
    clear i tes3_file data_file;
    %=====================================================================%
    
end