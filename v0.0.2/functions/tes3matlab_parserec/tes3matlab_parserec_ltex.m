function [ltex_name, ltex_intv, ...
          ltex_data, ltex_dele] = tes3matlab_parserec_ltex(raw_data)
      
    % Memory trackers
    n = numel(raw_data);
    m = 0;
    
    % Init
    ltex_dele = false;
    ltex_name = '';
    ltex_intv = uint32(0);
    ltex_data = '';
    
    % Parse data
    while(m<n)
        
        % Entry type and size
        entry_type = native2unicode( raw_data(m+(1:4))' );
        entry_size = bits2nums( raw_data(m+(5:8)) );
        m          = m + 8;
        
        % Process
        switch entry_type
            case 'NAME';    ltex_name = native2unicode( raw_data(m+(1:entry_size))' );
            case 'INTV';    ltex_intv = bits2nums( raw_data(m+(1:entry_size)) );
            case 'DATA';    ltex_data = native2unicode( raw_data(m+(1:entry_size))' );
            case 'DELE';    ltex_dele = true;
            otherwise;      error(['Unknown field in LTEX record: ' entry_type]);
        end
        m = m + entry_size;
        
    end
        
    % Clean strings by removing blanks and leading/trailing spaces
    ltex_name = lower(deblank( ltex_name ));
    ltex_data = lower(deblank( ltex_data ));
    
    % Remove leading directory path
    ii = strfind(ltex_data,'\');
    if(~isempty(ii));   ltex_data = ltex_data(ii(end)+1:end);   end
    
    % Remove file extension
    ii = strfind(ltex_data,'.');
    if(~isempty(ii));   ltex_data = ltex_data(1:(ii-1));  end
 
end