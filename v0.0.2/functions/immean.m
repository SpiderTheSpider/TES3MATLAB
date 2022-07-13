function y = immean(x)

    y = zeros(1,3,'uint8');
    for i=1:3
        z = double(x(:,:,i));
        y(i) = uint8(mean(z(:)));
    end

end