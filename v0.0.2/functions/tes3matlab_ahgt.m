function ahgt = tes3matlab_ahgt(vhgt,hoff)

    % Compute absolute vertex heights
    ahgt = zeros(size(vhgt),'int32');
	for ix=1:size(vhgt,1)
        ahgt(ix,1) = int32(hoff) + sum(int32(vhgt(1:ix,1)));
        for iy=2:size(vhgt,2)
        	ahgt(ix,iy) = ahgt(ix,iy-1) + int32(vhgt(ix,iy));
        end
    end
    
end