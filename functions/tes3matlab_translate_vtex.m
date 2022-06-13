function y = tes3matlab_translate_vtex(x)
% According to
% https://github.com/psi29a/tesannwyn/blob/master/tes3_export.c...
% 
% Create a sensible (X by Y) array from the strange MW format:
%  4x4 groups, bottom-left, bottom-right, top-left, top-right,
%  then those blocks re[ap]pear in a similar 4x4 pattern!

    % Init output and loop through quadrants
    y = zeros(16,16);
    for q=0:3
        
        % Quadrant
        switch q
            case 0; qx=1; qy=1; p=0;    % NOTE: MATLAB is 1-based, not 0-based like almost everyone else...
            case 1; qx=9; qy=1; p=32;
            case 2; qx=1; qy=9; p=128;
            case 3; qx=9; qy=9; p=160;
        end

        % Do
        for i=0:1
            p = p + (i*32);
            for j=0:1
            for k=0:3
            for l=0:3
            	y(p+1) = x(qy+k+(4*i),qx+l+(4*j));
                p = p + 1;
            end
            end
            end
        end
        clear i j k l;
      
    end
    
end