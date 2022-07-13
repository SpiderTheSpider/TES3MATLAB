function [cell_name, cell_rgnn, cell_flags, cell_x, cell_y, cell_nam5 ...
          ] = tes3matlab_parserec_cell(raw_data)

    % Memory tracker
    n = numel(raw_data);
    m = 0;
    
    % Init data
    cell_name  = '';                    % Name 
    cell_rgnn  = '';                    % Region
    cell_flags = [];                    % Flags
    cell_x     = int32(-inf);           % World grid x-position
    cell_y     = int32(-inf);           % World grid y-position
    cell_nam5  = zeros(1,4,'uint8');    % Map color (exterior only)
    
    % First two fields, NAME and DATA, are required
    entry_type = native2unicode( raw_data(m+(1:4))' );
    entry_size = bits2nums( raw_data(m+(5:8)) );
    m          = m + 8;
    if(strcmp(entry_type,'NAME'))
        cell_name = lower(deblank(native2unicode( raw_data(m+(1:entry_size))' )));
    else
        error(['First field of CELL should be NAME. Instead, found ' entry_type]);
    end
    m = m + entry_size;
    
    % DATA
    entry_type = native2unicode( raw_data(m+(1:4))' );
    entry_size = bits2nums( raw_data(m+(5:8)) );
    m          = m + 8;
    if(strcmp(entry_type,'DATA'))
        cell_flags = typecast( raw_data(m+(1:4) ) , 'uint8' );  
        cell_x     = typecast( raw_data(m+(5:8) ) , 'int32'  );
        cell_y     = typecast( raw_data(m+(9:12)) , 'int32'  );
    else
        error(['Second field of CELL should be DATA. Instead, found ' entry_type]);
    end
    m = m + entry_size;
    
    % Parse data
    while(m<n)
        
        % Entry type and size
        entry_type = native2unicode( raw_data(m+(1:4))' );
        entry_size = bits2nums( raw_data(m+(5:8)) );
        m          = m + 8;
        
        % Process
        switch entry_type
            %case 'NAME';    if(isempty(cell_name)); cell_name = lower(deblank( native2unicode( raw_data(m+(1:entry_size))' ) ));    end
            case 'RGNN';    cell_rgnn = lower(deblank( native2unicode( raw_data(m+(1:entry_size))' ) ));
            case 'NAM5';    cell_nam5 = raw_data(m+(1:entry_size));
            %case 'DATA'
            %    if(isempty(cell_flags))
            %        cell_flags = typecast( raw_data(m+(1:4) ) , 'uint8' );  
            %        cell_x     = typecast( raw_data(m+(5:8) ) , 'int32'  );
            %        cell_y     = typecast( raw_data(m+(9:12)) , 'int32'  );
            %    end
            otherwise;      %warning([' Skipping CELL field ' entry_type]); % Skip it
        end
        m = m + entry_size;
            
    end
    
    % Process data?
    
end