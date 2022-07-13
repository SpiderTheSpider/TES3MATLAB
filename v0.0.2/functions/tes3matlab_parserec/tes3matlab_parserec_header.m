function rec_data = tes3matlab_parserec_header(raw_data)

    % Init memory counters
    n = uint32(numel(raw_data));    % Total memory taken up by record
    m = uint32(0);                  % Current position in record memory
    
    % Init record data structures
    hedr = struct;
    mast = struct;
    
    % First 8 bits of header are empty
    m = m + 8;
            
    % First field should be HEDR
    field_name = native2unicode(raw_data(m+uint32(1:4))');
    if(strcmp(field_name,'HEDR')==false)
        error(['Expected first field of header to be HEDR, ' ...
               'instead, found ' field_name '!']);
    end
    m = m + 4;
    
    % HEDR size
    % NOTE: Instead of typecast, can do something like
    % b = a(1) + a(2)*256 + a(3)*(256^2) + a(4)*(256^3)
    hedr.size = typecast(raw_data(m+uint32(1:4)),'uint32');  
    m = m + 4;
            
    % File version (hmm, get odd results with this one...)
    hedr.vers = raw_data(m+uint32(1:4));
    m = m + 4;

    % File type (0=esm, 1=esp, 2=ess (script file))
    hedr.type = typecast(raw_data(m+uint32(1:4)),'uint32');
    m = m + 4;
            
    % File author name (can be an individual, team, or company)
    hedr.auth = native2unicode(raw_data(m+uint32(1:32)))';
    m = m + 32;
            
    % File description
    hedr.desc = native2unicode(raw_data(m+uint32(1:256)))';
    m = m + 256;
        
    % Number of records in file
    hedr.nrec = typecast(raw_data(m+uint32(1:4)),'uint32');
    m = m + 4;
            
    % If there is still more data, then they represent MAST fields
    n_mast = 0;
    while(m<n)
                
        % Should be a MAST field
        field_name = native2unicode(raw_data(m+uint32(1:4))');
        if(strcmp(field_name,'MAST')==false)
            error(['Expected a field with name MAST, instead ' ...
            'encountered ' field_name '!']);
        end
        m = m + 4;
        
        % A new master entry!
        n_mast = n_mast + 1;
        
        % Master file name length
        tmp = typecast(raw_data(m+uint32(1:4)),'uint32');
        m = m + 4;
        
        % Master file name
        mast(n_mast).name = native2unicode( raw_data(m+uint32(1:tmp)) )';
        m = m + tmp;
        
        % Remove case sensitivity and leading & trailing blanks
        mast(n_mast).name = lower(deblank( mast(n_mast).name ));
        
        % Skip 'DATA', data_size=8 entries
        m = m + 8;
        
        % Master file size
        mast(n_mast).size = typecast(raw_data(m+uint32(1:8)),'uint64');
        m = m + 8; 
     
    end %if(m<n)

    % Set outputs
    rec_data.hedr = hedr;
    rec_data.mast = mast;
    
end