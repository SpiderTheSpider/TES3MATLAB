function rec_data = tes3matlab_parserec(rec_type,raw_data)

    % Expect rec_data to be raw bits
    raw_data = uint8(raw_data);

    % Parse fields in record depending on the record type
    switch(rec_type)
        case 'TES3';    rec_data = tes3matlab_parserec_header(raw_data);
            
        otherwise;      error(['Unknown TES3 record type: ' rec_type]);
    end %switch(rec_type)
    
end