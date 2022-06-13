function med_color = median_color(input_image)
% Finds median color of image
    med_color = zeros(1,3);
    for i=1:3
        tmp = input_image(:,:,i);
        med_color(i) = median(tmp(:));
    end
end
