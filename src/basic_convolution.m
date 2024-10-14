function [filtered] = basic_convolution(image,kernel) 
    dimensions = size(image); 
    dimensions2 = size(kernel); 
 
    % define kernel center indices 
    kernelCenter_x = floor(dimensions2(1)/2); 
    kernelCenter_y = floor(dimensions2(2)/2); 
 
    image2 = zeros(dimensions(1),dimensions(2)); 
    for i = 1:dimensions(1) 
        for j = 1:dimensions(2) 
            for k = 1:dimensions2(1) 
                for l = 1:dimensions2(2) 
                    % New changes are added below 
                    ii = i+(k-kernelCenter_x); 
                    jj = j+(l-kernelCenter_y);  
                    if (ii >= 1 && ii <= dimensions(1) && jj >= 1 && jj <= dimensions(2)) 
                        %fprintf("i:%d j:%d ii:%d jj:%d k:%d l:%d\n",i,j,ii,jj,k,l);
                        image2(i,j) = image2(i,j) + image(ii,jj)* kernel(k,l); 
                    end 
                end 
            end 
        end   
    filtered = image2; 
    end 
end