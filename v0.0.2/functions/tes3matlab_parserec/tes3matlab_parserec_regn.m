function [regn_name, regn_fnam, regn_cnam, regn_snam, regn_weat, ...
          regn_bnam] = tes3matlab_parserec_regn(regn_data,n_weat)
      
    % Memory trackers
    n = numel(regn_data);
    m = 0;
      
    % Init outputs
    regn_name = '';                         % + ID
    regn_fnam = '';                         % + Full Name
    regn_cnam = zeros(1,4,'uint8');         % + Color 
    regn_snam = struct;                     % * Sound(s), name & chance
    regn_weat = zeros(1,n_weat,'uint8');	% + Weather chances       
    regn_bnam = '';                         % - Sleeping creature
      
    % Gather fields
    while(m<n)
        
        % Field entry name & size
        entry_name = native2unicode(regn_data(m+(1:4))');
        entry_size = bits2nums( regn_data(m+(5:8)) );
        m          = m + 8;
        
        % Field entry data
        entry_data = regn_data(m+(1:entry_size));
        m          = m + entry_size;
        
        % Parse field entry
        switch(entry_name)
        	case 'NAME';	regn_name = lower(deblank(native2unicode(entry_data')));
            case 'FNAM';	regn_fnam = lower(deblank(native2unicode(entry_data')));
            case 'CNAM';    regn_cnam = entry_data;
            case 'WEAT';    regn_weat(1:entry_size) = entry_data; 
            case 'BNAM';    regn_bnam = lower(deblank(native2unicode(entry_data')));
            case 'SNAM'   
                regn_snam(end+1).sound_name = lower(deblank(native2unicode(entry_data(1:32)')));
                regn_snam(end).sound_chance = entry_data(33);
            otherwise;      error(['  Unknown REGN sub-entry ' entry_name]);
        end %switch(entry_name)
        
    end %while(m<n)
      
end