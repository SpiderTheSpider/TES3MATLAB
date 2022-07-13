function file_data = tes3matlab_clean(file_data)
% NOTE: This should work for merged_data too, but it's probably better to
% clean before merging

    % Loop through files
    n_file = numel(file_data);
    for i=1:n_file

        %-----------------------------------------------------------------%
        % Remove empty region entries
        if(~isempty(fieldnames(file_data(i).regn)))
            j = 1;
            while(j<=numel(file_data(i).regn))
            if(isempty(file_data(i).regn(j).name))
                disp(['Cleaning empty REGN entry #' num2str(j) ' from ' file_data(i).file '...']);
                jj = 1:numel(file_data(i).regn);
                jj(j) = [];
                file_data(i).regn = file_data(i).regn(jj);
            else
                j = j + 1;
            end
            end
            clear j jj;
        end
        
        %-----------------------------------------------------------------%
        % Remove empty land entries, or fix bad ones
        if(~isempty(fieldnames(file_data(i).land)))
            j = 1;
            while(j<=numel(file_data(i).land))
            if(isempty(file_data(i).land(j).vhgt))
                disp(['Cleaning empty LAND entry #' num2str(j) ' from ' file_data(i).file '...']);
                jj = 1:numel(file_data(i).land);
                jj(j) = [];
                file_data(i).land = file_data(i).land(jj);
            else

% NOTE: THIS ASSUMED VTEX=0 LINKS TO _LAND_DEFAULT, WHICH IS NOT CORRECT
%                 % If vtex is empty, set it to 0's
%                 if(~isfield(file_data(i).land(j),'vtex'))
%                     disp(['Fixing missing VTEX in LAND entry #' num2str(j) ' from ' file_data(i).file '...']);
%                     file_data(i).land(j).vtex = zeros(16,16);
%                 elseif(isempty(file_data(i).land(j).vtex))
%                     disp(['Fixing empty VTEX in LAND entry #' num2str(j) ' from ' file_data(i).file '...']);
%                     file_data(i).land(j).vtex = zeros(16,16);
%                 end


                % Onto the next one
                j = j + 1;

            end
            end
            clear j jj;
        end
        
        %-----------------------------------------------------------------%
        % Remove unused ltex entries
%         if(~isempty(fieldnames(file_data(i).ltex)) & ~isempty([file_data(i).land.vtex])) 
%             
%             % Unique LTEX values
%             intv_defined = [0 [file_data(i).ltex.intv]+1];
%             intv_called = unique([file_data(i).land.vtex]);
%             ii = find(~ismember(intv_defined,intv_called)==1);
%         
%         end
        %-----------------------------------------------------------------%
        % Remove cell entries not referenced to a land
        
        
    end
    clear i;

end