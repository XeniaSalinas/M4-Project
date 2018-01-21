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
tmp_pad = padarray(L_image, [left_pad, left_pad], 'symmetric', 'pre');
L_image_pad = padarray(tmp_pad, [right_pad, right_pad], 'symmetric', 'post');

tmp_pad = padarray(R_image, [left_pad, left_pad], 'symmetric', 'pre');
R_image_pad = padarray(tmp_pad, [right_pad, right_pad], 'symmetric', 'post');
clear tmp_pad;
    

% Pre-allocate disparity map
[h, w] = size(L_image);
disp_map = zeros(h, w);

% Go over all pixels (i,j)
for i=left_pad+1:h+left_pad
    for j=left_pad+1:w+left_pad
        
        % Reference window values (left image)
        left_window_vals = L_image_pad(i-left_pad:i+right_pad,j-left_pad:j+right_pad);
        
        % Minimum and maximum disparities (left of origin)
        min_left_disparity_pos = max(j-max_disparity, left_pad+1);
        max_left_disparity_pos = j-min_disparity;
        % Minimum and maximum disparities (right of origin)
        min_right_disparity_pos = j+min_disparity;
        max_right_disparity_pos = min(j+max_disparity, w+left_pad);
        % Disparities interval
        if max_left_disparity_pos < j
            % Impossible configuration in left border, just use right
            % interval. Only happens when min_disparity is 0.
            interval = min_right_disparity_pos:max_right_disparity_pos;
        elseif min_right_disparity_pos > j
            % Impossible configuration in right border, just use left
            % interval. Only happens when min_disparity is 0.
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
        cost_vector = zeros(1, interval_idx);
        for idx=1:interval_idx
            % Column index of the moving window
            window_j = interval(idx);
            right_window_vals = R_image_pad(        ...
                i-left_pad:i+right_pad,             ...
                window_j-left_pad:window_j+right_pad  ...
            );
        
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
                gammap =17.5;
                T = 45;            
                left_vals = left_window_vals(:);
                right_vals = right_window_vals(:);
                center = ceil(window_size*window_size/2);
                left_center = left_vals(center);
                right_center = right_vals(center);             
                num = 0;
                den = 0;
                for k = 1:window_size*window_size
                    dist = sqrt((k-center)^2);
                    wleft = exp( - (abs(left_vals(k)-left_center)/gammac) - (dist/gammap));
                    wright = exp( - (abs(right_vals(k)-right_center)/gammac) - (dist/gammap));
                    c = min(abs(left_vals(k) - right_vals(k)),T);
                    num = num + wleft*wright*c;
                    den = den + wleft*wright;
                end
                bilateral_cost = num/den;
                cost_vector(idx) = bilateral_cost;
            end
        end
        
        % Pick disparity that minimizes cost / maximizes quality
        if strcmp(matching_cost, 'SSD') || strcmp(matching_cost, 'BILATERAL')
            [~, disp_pos] = min(cost_vector);
        elseif strcmp(matching_cost, 'NCC')
            [~, disp_pos] = max(cost_vector);
        end
        signed_disparity = interval(disp_pos) - j;
        disparity = abs(signed_disparity);
        
        % Assign disparity to disparity map
        disp_map(i-left_pad, j-left_pad) = disparity;
    end    
end


end

