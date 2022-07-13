function [land_dele, land_x, land_y, land_flag, land_vnml, land_hoff, ...
          land_vhgt, land_wnam, land_vclr, land_vtex ...
         ] = tes3matlab_parserec_land(raw_data,nv,nt,nw)

    % Memory tracker
    n = numel(raw_data);
    m = 0;

    % Init data
    land_dele = false;                  % Detected an issue, should delete?
    land_x    = int8(-inf);             % x-position on world (cell) map
    land_y    = int8(-inf);             % y-position on world (cell) map
    land_flag = zeros(1,4,'uint8');     % Flags, id valid fields present
    land_vnml = zeros(nv,nv,3,'int8');  % Vertex normals
    land_hoff = single(0);              % Entire cell height offset
    land_vhgt = zeros(nv,nv,'int8');    % Vertex heights
    land_wnam = zeros(nw,nw,'uint8');     % World minimap topography
    land_vclr = zeros(nv,nv,3,'uint8'); % Vertex colors
    land_vtex = zeros(nt,nt,'uint16');
    
    % Parse data
    while(m<n)

        % Entry type and size
        entry_type = native2unicode( raw_data(m+(1:4))' );
        entry_size = double(typecast( raw_data(m+(5:8)) , 'uint32' ));
        m          = m + 8;
        
        % Process
        switch entry_type
            case 'DATA';    land_flag = raw_data(m+(1:entry_size));
            case 'INTV'
                land_x = typecast( raw_data(m+(1:4)), 'int32' );
                land_y = typecast( raw_data(m+(5:8)), 'int32' );
            case 'VNML'
                mm = m;
                for ix=1:nv
                for iy=1:nv
                for iz=1:3
                    mm = mm + 1;
                    land_vnml(ix,iy,iz) = typecast( raw_data(mm), 'int8' );
                end
                end
                end
            case 'VHGT'
                mm = m;
                land_hoff = typecast( raw_data(mm+(1:4)), 'single' );
                mm = mm + 4;
                for ix=1:nv
                for iy=1:nv
                    mm = mm + 1;
                    land_vhgt(ix,iy) = typecast( raw_data(mm), 'int8' );
                end
                end
            case 'WNAM'
                mm = m;
                for ix=1:nw
                for iy=1:nw
                    mm = mm + 1;
                    land_wnam(ix,iy) = raw_data(mm);
                end
                end
            case 'VCLR'
                mm = m;
                for ix=1:nv
                for iy=1:nv
                for ic=1:3
                    mm = mm + 1;
                    land_vclr(ix,iy,ic) = raw_data(mm);
                end
                end
                end
            case 'VTEX'
                mm = m;
                for ix=1:nt
                for iy=1:nt
                    land_vtex(ix,iy) = typecast( raw_data(mm+(1:2)), 'uint16' );
                    mm = mm + 2;
                end
                end 
            otherwise;      %warning([' Skipping LAND field ' entry_type]); % Skip it
        end
        m = m + entry_size;
            
    end
    
    % Process data?
    if(~all(land_vtex(:)==0))
        land_vtex = tes3matlab_translate_vtex(land_vtex);
        land_vtex = land_vtex';
    end

end