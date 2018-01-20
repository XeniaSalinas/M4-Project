function disp_map = stereo_computation(L_image, R_image, min_disparity, ...
                                       max_disparity, window_size, ...
                                       matching_cost)
%STEREO_COMPUTATION Computes stereo disparity map from 2 images

% Assertions
if min_disparity < 0
    error('Minimum disparity should be 0 or more');
end

if window_size <= 0
    error('Window size should be bigger than 0');
end

if matching_cost ~= "SSD" && matching_cost ~= "NCC"
    error("Matching cost not recognized / implemented. Use one of: {'SSD', 'NCC'}");
end

% Convert images to greyscale
if length(size(L_image)) == 3
    L_image = sum(double(L_image), 3) / 3 / 255;
end
if length(size(R_image)) == 3
    R_image = sum(double(R_image), 3) / 3 / 255;
end

% Add padding to the image according to window_size
if mod(window_size, 2) == 1
    % Odd window size
    left_pad = floor(window_size / 2); right_pad = left_pad;
else
    % Even window size
    left_pad = window_size / 2 - 1; right_pad = window_size / 2;
end
tmp_pad = padarray(L_image, [left_pad, left_pad], 0, 'pre');
L_image_pad = padarray(tmp_pad, [right_pad, right_pad], 0, 'post');

tmp_pad = padarray(R_image, [left_pad, left_pad], 0, 'pre');
R_image_pad = padarray(tmp_pad, [right_pad, right_pad], 0, 'post');
    

% Pre-allocate disparity map
[h, w] = size(L_image);
disp_map = zeros(h, w);

% Go over all pixels (i,j)
for i=left_pad+1:h+left_pad
    for j=left_pad+1:w+left_pad
        
        % Compute window values, depends on the current pixel position
        % due to borders
%         window_idx = 
%         window_values_l = L_image_pad(i-
        
        % Compute vector of cost that depends on disparity
        
        % Pick disparity that minimizes cost / maximizes quality
        
    end    
end


end

