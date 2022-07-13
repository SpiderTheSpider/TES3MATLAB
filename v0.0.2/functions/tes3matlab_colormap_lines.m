function clrs = tes3matlab_colormap_lines(n)

    % Colormap: Lines
    clrs = zeros(3,379);
    clrs(:,001:064) = [zeros(1,64); zeros(1,64); 0:4:252];  % B
    clrs(:,065:127) = [zeros(1,63); 252:-4:4;    252:-4:4];    % C
    clrs(:,128:190) = [zeros(1,63); 4:4:252;     zeros(1,63)];  % G
    clrs(:,191:253) = [252:-4:4;    252:-4:4;    zeros(1,63)];    % Y
    clrs(:,254:316) = [4:4:252;     zeros(1,63); zeros(1,63)];  % R
    clrs(:,317:379) = [252:-4:4;    zeros(1,63); 252:-4:4];    % M

    % Convert
    clrs = (clrs)'./255;
    
    % Interpolate
    m    = size(clrs,1);
    clrs = interp1(1:m,clrs,1:(m/n):m,'linear');

end