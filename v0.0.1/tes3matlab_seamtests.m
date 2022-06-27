%function [ztest_results, ctest_results, ttest_results] = tes3matlab_seamtests(merged_data,map_data)
% log = tes3matlab_seamtests(map_data)
% Perform tests to look for seams, where data from adjacent cells aren't
% aligned, or where >2 land textures are present

    % List of landscape cells
    xy = [merged_data.land.xy];
    x  = unique(xy(1,:));   nx = numel(x);  x0 = min(x);
    y  = unique(xy(2,:));   ny = numel(y);  y0 = min(y);
    n  = size(xy,2);
    
    % Loop through land entries
    ztest_log = zeros(0,66);     ztest_img = zeros(nx*8,ny*8,'uint8');
    ctest_log = zeros(0,66);     ctest_img = zeros(nx*8,ny*8,'uint8');
    ttest_log = zeros(0,18);     ttst_img  = zeros(nx*16,ny*16,'uint8');
    
    % Pairs that have been tested
    pair_list_z = 0;
    pair_list_c = 0;
    
    % Loop
    for i=1:n

    	% This cell's location on the world map
        ix = xy(1,i)-x0+1;     %ix = (ix-1)*8+1 : ix*8;
        iy = xy(2,i)-y0+1;     %iy = (iy-1)*8+1 : iy*8;

        % Look for cells to left, right, top, and bottom
        i_lft = find(xy(1,:)==xy(1,i)-1 & xy(2,:)==xy(2,i)  );  if(isempty(i_lft));	i_lft = 0;  end  
        i_rgt = find(xy(1,:)==xy(1,i)+1 & xy(2,:)==xy(2,i)  );  if(isempty(i_rgt));	i_rgt = 0;  end 
        i_btm = find(xy(1,:)==xy(1,i)   & xy(2,:)==xy(2,i)-1);  if(isempty(i_btm));	i_btm = 0;  end 
        i_top = find(xy(1,:)==xy(1,i)   & xy(2,:)==xy(2,i)+1);  if(isempty(i_top));	i_top = 0;  end 
        
        %-----------------------------------------------------------------%
        % HEIGHTS
        if(~isempty(merged_data.land(i).ahgt))
    
            % This cell's span in the log image
            ii = (ix-1)*8+1 : ix*8;
            jj = (iy-1)*8+1 : iy*8;

            % This cell's absolute vertex heights
            z0 = merged_data.land(i).ahgt';

            % Gather height datums if they exist
            % Else, assume it's an undefined cell (all heights are -256)
            if(i_lft==0);	z_lft = -256.*ones(65,65,'int32');	else;   z_lft = merged_data.land(i_lft).ahgt';	end
            if(i_rgt==0);   z_rgt = -256.*ones(65,65,'int32');  else;   z_rgt = merged_data.land(i_rgt).ahgt';	end
            if(i_btm==0);   z_btm = -256.*ones(65,65,'int32');  else;   z_btm = merged_data.land(i_btm).ahgt';	end
            if(i_top==0);   z_top = -256.*ones(65,65,'int32');  else;   z_top = merged_data.land(i_top).ahgt';	end

            % NOTE: DUE TO HOW IMAGES ARE STORED, FULL IMAGE IS ACTUALLY
            %       lft
            % btm   000   top
            %       rgt   
            % THEN THAT TRANSPOSED

            % To visually confirm this, look at the following:
            %zz = [z0(:) z_lft(:) z_rgt(:) z_btm(:) z_top(:)];
            %zr = [min(zz(:)) max(zz(:))];
            %subplot(3,3,2); imshow(z_lft,[]);   caxis(zr);
            %subplot(3,3,4); imshow(z_btm,[]);   caxis(zr);
            %subplot(3,3,5); imshow(z0,[]);      caxis(zr);
            %subplot(3,3,6); imshow(z_top,[]);   caxis(zr);
            %subplot(3,3,8); imshow(z_rgt,[]);   caxis(zr);

            % Check seam between lft and 0
            check1 = z_lft(end,1:64);
            check2 = z0(1,1:64);
            if(~all(check1==check2))
                if(i_lft==0); a=-1; else; a=i_lft; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_z,b) & ~ismember(pair_list_z,c))
                    ztest_log(end+1,:) = [i, a, abs(check1(:)-check2(:))'];
                    pair_list_z(end+1) = b;
                    pair_list_z(end+1) = c;
                    ztest_img(ii(1),jj,:) = 1;
                end
            end

            % Check seam between btm and 0
            check1 = z_btm(1:64,end);
            check2 = z0(1:64,1);
            if(~all(check1==check2))
                if(i_btm==0); a=-2; else; a=i_btm; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_z,b) & ~ismember(pair_list_z,c))
                    ztest_log(end+1,:) = [i, a, abs(check1(:)-check2(:))'];
                    pair_list_z(end+1) = b;
                    pair_list_z(end+1) = c;
                    ztest_img(ii,jj(1),:) = 1;
                end
            end

            % Check seam between rgt and 0
            check1 = z_rgt(1,1:64);
            check2 = z0(end,1:64);
            if(~all(check1==check2))
                if(i_rgt==0); a=-3; else; a=i_rgt; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_z,b) & ~ismember(pair_list_z,c))
                    ztest_log(end+1,:) = [i, a, abs(check1(:)-check2(:))'];
                    pair_list_z(end+1) = b;
                    pair_list_z(end+1) = c;
                    ztest_img(ii(end),jj,:) = 1;
                end
            end

            % Check seam between top and 0
            check1 = z_top(1:64,1);
            check2 = z0(1:64,end);
            if(~all(check1==check2))
                if(i_top==0); a=-4; else; a=i_top; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_z,b) & ~ismember(pair_list_z,c))
                    ztest_log(end+1,:) = [i, a, abs(check1(:)-check2(:))'];
                    pair_list_z(end+1) = b;
                    pair_list_z(end+1) = c;
                    ztest_img(ii,jj(end),:) = 1;
                end
            end
        
        end
        %-----------------------------------------------------------------%
        % COLORS
        if(~isempty(merged_data.land(i).vclr))
            
            % This cell's span in the log image
            ii = (ix-1)*8+1 : ix*8;
            jj = (iy-1)*8+1 : iy*8;    
            
            % This cell's vertex colors
            c0 = permute(merged_data.land(i).vclr,[2 1 3]);
    
            % Gather vertex color datums if they exist
            c_lft = []; if(i_lft~=0);   if(~isempty(merged_data.land(i_lft).vclr));	c_lft = permute(merged_data.land(i_lft).vclr,[2 1 3]);	end;    end
            c_rgt = []; if(i_rgt~=0);   if(~isempty(merged_data.land(i_rgt).vclr));	c_rgt = permute(merged_data.land(i_rgt).vclr,[2 1 3]);	end;    end
            c_btm = []; if(i_btm~=0);   if(~isempty(merged_data.land(i_btm).vclr));	c_btm = permute(merged_data.land(i_btm).vclr,[2 1 3]);	end;    end
            c_top = []; if(i_top~=0);   if(~isempty(merged_data.land(i_top).vclr));	c_top = permute(merged_data.land(i_top).vclr,[2 1 3]);	end;    end

            % Check seam between lft and 0
            if(~isempty(c_lft))
                if(i_lft==0); a=-1; else; a=i_lft; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_c,b) & ~ismember(pair_list_c,c))    
                    check1 = squeeze(c_lft(end,1:64,:));
                	check2 = squeeze(c0(1,1:64,:));
                    test   = [ check1(:,1)~=check2(:,1) | check1(:,2)~=check2(:,2) | check1(:,3)~=check2(:,3) ];
                    if(~all(test==0))
                        ctest_log(end+1,:) = [i,a,test(:)'];
                        pair_list_c(end+1) = b;
                        pair_list_c(end+1) = c;
                        ctest_img(ii(1),jj,:) = 1;
                    end
                end
            end
            
            % Check seam between btm and 0
            if(~isempty(c_btm))
            	if(i_btm==0); a=-2; else; a=i_btm; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_c,b) & ~ismember(pair_list_c,c))        
                    check1 = squeeze(c_btm(1:64,end,:));
                    check2 = squeeze(c0(1:64,1,:));
                    test   = [ check1(:,1)~=check2(:,1) | check1(:,2)~=check2(:,2) | check1(:,3)~=check2(:,3) ];
                    if(~all(test==0))
                        ctest_log(end+1,:) = [i,a,test(:)'];
                        pair_list_c(end+1) = b;
                        pair_list_c(end+1) = c;
                        ctest_img(ii,jj(1),:) = 1;
                    end
                end
            end
            
            % Check seam between rgt and 0
            if(~isempty(c_rgt))
                if(i_rgt==0); a=-3; else; a=i_rgt; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_c,b) & ~ismember(pair_list_c,c))    
                    check1 = squeeze(c_rgt(1,1:64,:));
                    check2 = squeeze(c0(end,1:64,:));
                    test   = [ check1(:,1)~=check2(:,1) | check1(:,2)~=check2(:,2) | check1(:,3)~=check2(:,3) ];
                    if(~all(test==0))
                        ctest_log(end+1,:) = [i,a,test(:)'];
                        pair_list_c(end+1) = b;
                        pair_list_c(end+1) = c;
                        ctest_img(ii(end),jj,:) = 1;
                    end
                end
            end
            
            % Check seam between top and 0
            if(~isempty(c_top))
                if(i_top==0); a=-4; else; a=i_top; end
                b = i + sqrt(-1)*a;
                c = sqrt(-1)*i + a;
                if(~ismember(pair_list_c,b) & ~ismember(pair_list_c,c))    
                    check1 = squeeze(c_top(1:64,1,:));
                    check2 = squeeze(c0(1:64,end,:));
                    test   = [ check1(:,1)~=check2(:,1) | check1(:,2)~=check2(:,2) | check1(:,3)~=check2(:,3) ];
                    if(~all(test==0))
                        ctest_log(end+1,:) = [i,a,test(:)'];
                        pair_list_c(end+1) = b;
                        pair_list_c(end+1) = c;
                        ctest_img(ii,jj(end),:) = 1;
                    end
                end
            end

        end
        %-----------------------------------------------------------------%
        % TEXTURES
        %if(~isempty(merged_data.land(i).vclr))   
        %    % ACTUALLY, JUST USE FULL LTEX IMAGE 
        %end
            
    end
    clear i;
    
%     % TEXTURES TEST [NOW DONE IN ANOTHER CODE]
%     ttest_img = zeros(size(map_data.ltex,1)-1,size(map_data.ltex,2)-1,'uint8');
%     for i=1:size(map_data.ltex,1)-1
%     for j=1:size(map_data.ltex,2)-1
%        
%         % Surrounding 4 vertex colors
%         tmp = map_data.ltex(i:(i+1),j:(j+1));
%         ttest_img(i,j) = numel(unique(tmp(:)));
%         
%     end
%     end
    
    % Permute images
    ztest_img = permute(ztest_img(:,end:-1:1,:),[2 1 3]);
    ctest_img = permute(ctest_img(:,end:-1:1,:),[2 1 3]);
% ttest_img = permute(ttest_img(:,end:-1:1,:),[2 1 3]);
    
    % Write logs as text files
    imwrite(uint8(ztest_img.*255),[opts.dir_xprt '\seamtest_vhgt.png'],'png')
    imwrite(uint8(ctest_img.*255),[opts.dir_xprt '\seamtest_vtex.png'],'png')
    
        
%end