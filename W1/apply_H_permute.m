function I2 = apply_H(I, H)
%APPLY_H Transform the input image with the given homography
%% Compute the size of the transformed image

I = permute(I,[2 1 3]);

[n,m,c] = size(I);
% Homogeneous coordinates of the corners
lu_p = double([0, 0, 1]');  % Left-up
ru_p = double([0, m, 1]');  % Right-up
ld_p = double([n, 0, 1]');  % Left-down
rd_p = double([n, m, 1]');  % Right down

% Homogeneous coordinates of the transformed image's corners
lu2_p = H * lu_p;  % Left-up 
ru2_p = H * ru_p;  % Right-up
ld2_p = H * ld_p;  % Left-down
rd2_p = H * rd_p;  % Right down

homogeneous_corners = [lu2_p'; ru2_p'; ld2_p'; rd2_p'];
    
% Cartesian coordinates of the transformed image's corners
cartesian_corners = [homogeneous_corners(:, 1)./homogeneous_corners(:, 3) ...
                        homogeneous_corners(:, 2)./homogeneous_corners(:, 3)];

% top left and bottom right new corners
lu_c = min(cartesian_corners);
rd_c = max(cartesian_corners);

% Rows (m2) and columns (n2) of the transformed image
n2 = round(rd_c(1) - lu_c(1));
m2 = round(rd_c(2) - lu_c(2));

%% Compute the homography transformation
I2 = nan(n2, m2, c);

% Variables to store the output coordinates
indices = zeros(n2,m2);
x2_aux_1(:,:,1) = nan(n2,m2);
x2_aux_2(:,:,2) = nan(n2,m2);
        
for i = 1:n2
    for j = 1:m2
        x2_e = [i; j] + lu_c';
        x2_p = H\[x2_e;1]; %A\b for inv(a)*b  
        x_e = [x2_p(1)/x2_p(3), x2_p(2)/x2_p(3)];
        % If the transformed point is inside the source image, save the
        % point coordinates
        if (x_e(1) > 1 && x_e(2) > 1 && x_e(1) < n && x_e(2) < m)
            indices(i,j) = 1;
            x2_aux_1(i,j) = x_e(1);
            x2_aux_2(i,j) = x_e(2);
        end
    end
end

% Initialize interpolated image with NaN
I_interp = nan(numel(I2(:,:,1)),3);

% Coordinates of the source image
[X,Y]=meshgrid(1:m,1:n);

% Interpolate each channel
I_interp(indices>0,1) = interp2(X,Y, double(I(:,:,1)),x2_aux_2(indices>0),x2_aux_1(indices>0));
I_interp(indices>0,2) = interp2(X,Y, double(I(:,:,2)),x2_aux_2(indices>0),x2_aux_1(indices>0));
I_interp(indices>0,3) = interp2(X,Y, double(I(:,:,3)),x2_aux_2(indices>0),x2_aux_1(indices>0));

% Reshape the interpolated image
I2 = reshape(I_interp,n2,m2,3);
% Convert NaN to zeros
I2(isnan(I2))=0;

I2 = permute(I2,[2 1 3]);

end

