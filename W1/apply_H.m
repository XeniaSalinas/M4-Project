function I2 = apply_H(I, H)
%APPLY_H Transform the input image with the given homography

%% Compute the size of the transformed image
[n,m,c] = size(I);

% Cartesian coordinates of the corners
xmin = 1;
xmax = m;
ymin = 1;
ymax = n;

lu_c = [xmin; ymin];  % Left-up
ru_c = [xmax; ymin];  % Right-up
ld_c = [xmin; ymax];  % Left-down
rd_c = [xmax; ymax];  % Right down

% Homogeneous coordinates of the corners
lu_p = [lu_c;1];  
ru_p = [ru_c;1];  
ld_p = [ld_c;1];  
rd_p = [rd_c;1];  

% Homogeneous coordinates of the transformed image's corners
lu2_p = H * lu_p;  % Left-up 
ru2_p = H * ru_p;  % Right-up
ld2_p = H * ld_p;  % Left-down
rd2_p = H * rd_p;  % Right down

% Cartesian coordinates of the transformed image's corners
lu2_c = [lu2_p(1)/lu2_p(3), lu2_p(2)/lu2_p(3)];
ru2_c = [ru2_p(1)/ru2_p(3), ru2_p(2)/ru2_p(3)];
ld2_c = [ld2_p(1)/ld2_p(3), ld2_p(2)/ld2_p(3)];
rd2_c = [rd2_p(1)/rd2_p(3), rd2_p(2)/rd2_p(3)];

% % Size of the transformed image
% m2 = round(abs(max([lu2_c(1), ru2_c(1), ld2_c(1), rd2_c(1)]) - min([lu2_c(1), ru2_c(1), ld2_c(1), rd2_c(1)])));
% n2 = round(abs(max([lu2_c(2), ru2_c(2), ld2_c(2), rd2_c(2)]) - min([lu2_c(2), ru2_c(2), ld2_c(2), rd2_c(2)])));
% 
% % Compute the position rectification vector
% xoffset = floor(min([lu2_c(1), ru2_c(1), ld2_c(1), rd2_c(1)]));
% yoffset = floor(min([lu2_c(2), ru2_c(2), ld2_c(2), rd2_c(2)]));
             
% top left and bottom right new corners
lu_c = min([lu2_c; ru2_c; ld2_c; rd2_c]);
rd_c = max([lu2_c; ru2_c; ld2_c; rd2_c]);

% Rows (m2) and columns (n2) of the transformed image
m2 = round(rd_c(1) - lu_c(1));
n2 = round(rd_c(2) - lu_c(2));

%% Compute the homography transformation
I2 = zeros(n2,m2,c);
% Variables to store the output coordinates
indices = zeros(n2,m2);
x2_aux_1 = zeros(n2,m2);
x2_aux_2 = zeros(n2,m2);

for i = 1:n2
    for j = 1:m2
        % Alias for better understanding
        x_idx = j;
        y_idx = i;
        
        % Cartesian coordinate of the transformed image to be obtained in
        % the original image
        %x2_e = [x_idx; y_idx] + [xoffset; yoffset];
        x2_e = [x_idx; y_idx] + lu_c';
        x2_p = H \ [x2_e; 1]; % same as inv(H)*b, but more accurate and efficient
        x_e = [x2_p(1)/x2_p(3), x2_p(2)/x2_p(3)];
        % If the transformed point is inside the source image, save the
        % point coordinates
        if (x_e(1) >= 1 && x_e(2) >= 1 && x_e(1) <= m && x_e(2) <= n)
            indices(i,j) = 1;
            x2_aux_1(i,j) = x_e(1);
            x2_aux_2(i,j) = x_e(2);
        end
    end
end

% Initialize interpolated image
I_interp = zeros(numel(I2(:,:,1)),3);

% Interpolate each channel
I_interp(indices>0,1) = interp2(double(I(:,:,1)),x2_aux_1(indices>0),x2_aux_2(indices>0));
I_interp(indices>0,2) = interp2(double(I(:,:,2)),x2_aux_1(indices>0),x2_aux_2(indices>0));
I_interp(indices>0,3) = interp2(double(I(:,:,3)),x2_aux_1(indices>0),x2_aux_2(indices>0));

% Reshape the interpolated image
I2 = reshape(I_interp,n2,m2,3);

end

