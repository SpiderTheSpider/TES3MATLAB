function file_info = tes3matlab_getfileinfo(file_path)
% file_info = tes3matlab_getfileinfo(file_path)
% Reads in and returns information about the TES3 file linked to by
% file_path, including:
%  header data (file information, number of records, master files (if any))
%  defined record types
%  number of each record type present
% 
% by Spider the spider
% on 07/01/2022

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DEFINITIONS
    
    % Variable definitions
    % NOTE: This step is not necessary, but it acts as a sort of glossary
    % for this code.  It may also help when converting this code into other
    % languages which define variables in a preamble (e.g. C languages).
    n_types     = uint8(42);        % # possible record types in TES3 file
    rec_types   = cell(n_types,1);  % Cell array of strings of record types
    fid         = [];               % File identifier
    header_name = char(zeros(4,1)); % Should be TES3
    header_size = uint32(0);        % Size of header record in bytes
    field_name  = char(zeros(4,1)); % Current record field name
    mem_on      = uint32(0);        % For keeping track of position in file
    nul         = [];               % For skipping data
    tmp         = [];               % For temporarily-used data
    
    % File_info data structure
    file_info = struct;             % File info
     file_info.hedr = struct;       % Header record HEDR field
      file_info.hedr.vers = [];         % File version (?)
      file_info.hedr.type = uint32(0);  % 0=esm, 1=esp, 2=ess?
      file_info.hedr.auth = char(zeros(32,1));  % File author/company
      file_info.hedr.desc = char(zeros(256,1)); % File description
     file_info.mast = struct;       % Header record MAST field(s)
      file_info.mast.name = '';     % Name of master file
      file_info.mast.size = [];     % Size of master file (in bytes)
     file_info.rec = struct;        % File records
      file_info.rec.type  = '';     % Record type (string)
      file_info.rec.size  = [];     % Record size (in bytes)
      file_info.rec.flag  = [];     % Record flags (4 bits)
      file_info.rec.pif   = [];     % Record position in file
     file_info.rec_types = {};      % Unique record types in this file
     file_info.rec_amnts = [];      % Number of each unique record type
     
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INIT

    % Possible record types
    % NOTE: MAST records are technically sub-records of the header, but 
    % aren't included in its SIZE field, and neither are they counted as
    % individual records, so they must be handled carefully...
    rec_types = tes3matlab_recdef;
    n_types   = numel(rec_types);

	% Open the file for reading
	fid = fopen(file_path,'rb');
    if(fid==-1); error(['Could not open ' file_path ' for reading!']);  end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARSE HEADER
    
	% File signature, also used as name of the Header record
	header_name = native2unicode(fread(fid,4,'ubit8')');
    if(strcmp(header_name,'TES3')~=true);
        error(['Expected ' input_file ' file header to be TES3, ' ...
               'instead it was ' header_name '!']);
    end
    
	% Get (remainder of) header record size
    % NOTE: includes size of HEDR + MAST fields
	header_size = fread(fid,1,'uint32')+8;  % Since there's 8 empty bits
    header_data = fread(fid,header_size,'ubit8');
    tmp = tes3matlab_parserec('TES3',header_data);
    
    % Save header info to output structure
    file_info.hedr = tmp.hedr;
    file_info.mast = tmp.mast;

    % Debug log
    hedr_nrec = file_info.hedr.nrec;
    disp(['Gathering info from ' num2str(hedr_nrec) ' records in file ' file_path '...']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MAIN LOOP
    
    % Init record containers
	rec_type = char(zeros(hedr_nrec,4));
    rec_pif  = zeros(hedr_nrec,1,'uint32'); % Position in file
    rec_size = zeros(hedr_nrec,1,'uint32');
    rec_flag = zeros(hedr_nrec,4,'uint8');
    %rec_data = cell(hedr_nrec,1);          % Nah, read data in later...
    
    % Loop through entries
    n_rec = 0;
    while(~feof(fid))

        % Record type
        this_type = native2unicode( fread(fid,4,'ubit8')' );
        
        % Error check: Encountered EOF
        if(feof(fid) && n_rec<hedr_nrec)
            error(['Unexpected EOF at record #' num2str(i_rec) ...
                   ', expected ' num2str(hedr_nrec) ' records!']);
        end

        %=================================================================%
        % Continue if record type is valid
        % If it's a MAST entry, parse it separately
        if(numel(this_type)==4)  % Get an empty this_type at EOF
        if(ismember(this_type,rec_types)==true)
            
            % A new record!
            n_rec = n_rec + 1;
            
            % Record info
            rec_pif(n_rec)    = ftell(fid)-4;           % Position in file
            rec_type(n_rec,:) = this_type;              % Record type
            rec_size(n_rec)   = fread(fid,1,'uint32');  % Record size
            fseek(fid,4,0);                             % Unused data
            rec_flag(n_rec,:) = fread(fid,4,'ubit8');   % Record flags
            fseek(fid,rec_size(n_rec),0);               % Skip remainder
            
        else
            
            % Error found: Unknown record type
            error(['Unknown record type at record #' num2str(n_rec) ...
                  rec_type '!']);
              
        end %if(ismember(rec_type,rec_types)==true)
        end %if(numel(this_type)>4)
        %=================================================================%
        
    end %while(~feof(fid))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EPILOGUE
    
    % Error check: Found all records?
    if(n_rec ~= hedr_nrec)
        warning(['File specified there are ' num2str(hedr_nrec) ' records,' ...
                'but only ' num2str(n_rec) ' were found!']);
    end
    
    % Compute number of each record type found
    n_type = zeros(n_types,1);
    for i=1:n_rec
        ii = find(strcmp(rec_types,rec_type(i,:))==true);
        n_type(ii) = n_type(ii) + 1;
    end
    
    % Write results to output structure
    for i=1:n_rec
        file_info.rec(i).type = rec_type(i,:);
        file_info.rec(i).size = rec_size(i);
        file_info.rec(i).flag = rec_flag(i,:);
        file_info.rec(i).pif  = rec_pif(i);
    end
    
    % Record types & number of each record type
    ii = find(n_type~=0);
    file_info.rec_types = rec_types(ii);
    file_info.rec_amnts = n_type(ii);

end