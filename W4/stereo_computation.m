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

if strcmp(matching_cost, 'SSD') && strcmp(matching_cost, 'NCC') && strcmp(matching_cost, 'BILATERAL')
    error('Matching cost not recognized / implemented. Use one of: {''SSD, NCC''}');
end
if strcmp(matching_cost, 'BILATERAL')
    L_image = rgb2lab(L_image);
    R_image = rgb2lab(R_image);
else
    % Convert images to greyscale
    if length(size(L_image)) == 3
        L_image = sum(double(L_image), 3) / 3 / 255;
    end
    if length(size(R_image)) == 3
        R_image = sum(double(R_image), 3) / 3 / 255;
    end
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
clear tmp_pad;
    

% Pre-allocate disparity map
[h, w, c] = size(L_image);
disp_map = zeros(h, w, c);

if strcmp(matching_cost, 'BILATERAL')
    positions = 1:window_size;
    center = ceil(length(positions)/2);
    centered_positions = positions - center;
    % We normalize the distance in pixels to be in a range between [-1,1]
    % for a given dimension.
    centered_positions_norm = centered_positions / (window_size * 2);
    column_distances = repmat(centered_positions_norm, window_size, 1);
    row_distances = transpose(column_distances);
    dist = zeros(window_size, window_size, 3);
    dist(:,:,1) = sqrt(column_distances.^2 + row_distances.^2);
    dist(:,:,2) = sqrt(column_distances.^2 + row_distances.^2);
    dist(:,:,3) = sqrt(column_distances.^2 + row_distances.^2);
end

% Go over all pixels (i,j)
for i=left_pad+1:h+left_pad
    for j=left_pad+1:w+left_pad
        
        % Reference window values (left image)
        left_window_vals = L_image_pad(i-left_pad:i+right_pad,j-left_pad:j+right_pad,:);
        
        % Minimum and maximum disparities (left of origin)
        min_left_disparity_pos = max(j-max_disparity, left_pad+1);
        max_left_disparity_pos = j-min_disparity;
        % Minimum and maximum disparities (right of origin)
        min_right_disparity_pos = j+min_disparity;
        max_right_disparity_pos = min(j+max_disparity, w+left_pad);
        % Disparities interval
        if max_left_disparity_pos < min_left_disparity_pos
            % Impossible configuration in left border, just use right
            % interval. Only happens when min_disparity is different than 0.
            interval = min_right_disparity_pos:max_right_disparity_pos;
        elseif min_right_disparity_pos > max_right_disparity_pos
            % Impossible configuration in right border, just use left
            % interval. Only happens when min_disparity is different than 0.
            interval = min_left_disparity_pos:max_left_disparity_pos;
        elseif max_left_disparity_pos == min_right_disparity_pos
            % Central pixel positions (not close to borders) when 
            % min_disparity is 0.
            left_interval = min_left_disparity_pos:max_left_disparity_pos - 1;
            right_interval = min_right_disparity_pos:max_right_disparity_pos;
            interval = [left_interval right_interval];
        else
            % Central pixel positions (not close to borders) when 
            % min_disparity is different than 0.
            left_interval = min_left_disparity_pos:max_left_disparity_pos - 1;
            right_interval = min_right_disparity_pos:max_right_disparity_pos;
            interval = [left_interval right_interval];
        end
        
        % Pre-allocate cost vector from which the disparity will extracted
        % (as a minimum if it uses a cost function or as a maximum if it
        % uses a quality criterion)
        interval_idx = length(interval);
        cost_vector = zeros(interval_idx, c);
        for idx=1:interval_idx
            % Column index of the moving window
            window_j = interval(idx);
            right_window_vals = R_image_pad(i-left_pad:i+right_pad,window_j-left_pad:window_j+right_pad,:);
        
            % Compute cost between left (reference) window and right
            % (sliding) window
            if strcmp(matching_cost, 'SSD')
                % Assume uniform distribution of weights w(p,q)
                left_vals = left_window_vals(:);
                right_vals = right_window_vals(:);
                ssd = sum((left_vals - right_vals).^2) / (window_size * window_size);
                cost_vector(idx) = ssd;
            elseif strcmp(matching_cost, 'NCC')
                % Assume uniform distribution of weights w(p,q)
                left_vals = left_window_vals(:);
                right_vals = right_window_vals(:);
                m_left = mean(left_vals);
                m_right = mean(right_vals);
                std_left = std(left_vals);
                std_right = std(right_vals);
                ncc_num = sum((left_vals - m_left).*(right_vals - m_right)) / (window_size * window_size);
                ncc = ncc_num / (std_left * std_right);
                cost_vector(idx) = ncc;
            elseif strcmp(matching_cost, 'BILATERAL')
                gammac = 5;
                gammap = window_size/2;
                left_center = left_window_vals(floor(window_size/2)+1, floor(window_size/2)+1,:);
                right_center = right_window_vals(floor(window_size/2)+1, floor(window_size/2)+1,:);             
                bw_left = exp( - (abs(left_window_vals-left_center)/gammac) - (dist/gammap));
                bw_right = exp( - (abs(right_window_vals-right_center)/gammac) - (dist/gammap));
                bw = bw_left.*bw_right;
                cost_vector(idx, 1) = (sum(sum(bw(:,:,1).*min(abs(left_window_vals(:,:,1)-right_window_vals(:,:,1)),max_disparity))))/(sum(sum(bw(:,:,1))));
                cost_vector(idx, 2) = (sum(sum(bw(:,:,2).*min(abs(left_window_vals(:,:,2)-right_window_vals(:,:,2)),max_disparity))))/(sum(sum(bw(:,:,2))));
                cost_vector(idx, 3) = (sum(sum(bw(:,:,3).*min(abs(left_window_vals(:,:,3)-right_window_vals(:,:,3)),max_disparity))))/(sum(sum(bw(:,:,3))));
            end
        end        
        % Pick disparity that minimizes cost / maximizes quality
        if strcmp(matching_cost, 'SSD') 
            [~, disp_pos] = min(cost_vector);
            signed_disparity = interval(disp_pos) - j;
            disparity = abs(signed_disparity);

            % Assign disparity to disparity map
            disp_map(i-left_pad, j-left_pad) = disparity;
        elseif strcmp(matching_cost, 'BILATERAL')
            [~, disp_pos1] = min(cost_vector(:,1));
            [~, disp_pos2] = min(cost_vector(:,2));
            [~, disp_pos3] = min(cost_vector(:,3));
            signed_disparity1 = interval(disp_pos1) - j;
            disparity1 = abs(signed_disparity1);
            signed_disparity2 = interval(disp_pos2) - j;
            disparity2 = abs(signed_disparity2);
            signed_disparity3 = interval(disp_pos3) - j;
            disparity3 = abs(signed_disparity3);
            disp_map(i-left_pad, j-left_pad, 1) = disparity1;
            disp_map(i-left_pad, j-left_pad, 2) = disparity2;
            disp_map(i-left_pad, j-left_pad, 3) = disparity3;
        elseif strcmp(matching_cost, 'NCC')
            [~, disp_pos] = max(cost_vector);
            signed_disparity = interval(disp_pos) - j;
            disparity = abs(signed_disparity);

            % Assign disparity to disparity map
            disp_map(i-left_pad, j-left_pad) = disparity;
        end
        
    end    
end


end

