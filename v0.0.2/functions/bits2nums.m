function y = bits2nums(b)

    % Translate every 4 bits into an integer
    n = numel(b);
    y = zeros(1,n/4);
    for i=1:n/4
       ii = (i-1)*4+1 : i*4;
       y(i) = b(ii(1)) + b(ii(2))*256 + b(ii(3))*(256^2) + b(ii(4)).*(256^3);
    end

end