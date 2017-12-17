function I2 = apply_H(I, H)
%APPLY_H Transform the input image with the given homography

%% Compute the size of the transformed image
[n,m,c] = size(I);
% Homogeneous coordinates of the corners
lu_p = [1;1;1];  % Left-up
ru_p = [1;m;1];  % Right-up
ld_p = [n;1;1];  % Left-down
rd_p = [n;m;1];  % Right down
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
% Rows (m2) and columns (n2) of the transformed image
n2 = ceil(abs(max([lu2_c(1), ru2_c(1), ld2_c(1), rd2_c(1)]) - min([lu2_c(1), ru2_c(1), ld2_c(1), rd2_c(1)]))) + 1;
m2 = ceil(abs(max([lu2_c(2), ru2_c(2), ld2_c(2), rd2_c(2)]) - min([lu2_c(2), ru2_c(2), ld2_c(2), rd2_c(2)]))) + 1;

%% Compute the position rectification vector
origin_n = min([lu2_c(1), ru2_c(1), ld2_c(1), rd2_c(1)]);
origin_m = min([lu2_c(2), ru2_c(2), ld2_c(2), rd2_c(2)]);

%% Compute the homography transformation
% Query points in X and Y axis for interpolation. By default we set the
% value to -1 so that it will return NaN (for undefined points, outside of
% the transformed image), and then we convert NaN to zeroes. 
Xq = -1 * ones(n2, m2);
Yq = -1 * ones(n2, m2);
I2 = zeros(n2, m2, c);

% Variables to store the output coordinates
indices = zeros(n2,m2);
x2_aux_1(:,:,1) = zeros(n2,m2);
x2_aux_2(:,:,2) = zeros(n2,m2);

for i = 1:n2
    for j = 1:m2
        x2_e = [i; j] + [origin_n; origin_m];
<<<<<<< HEAD
        x2_p = H\[x2_e;1]; %A\b for inv(a)*b  
        x_e = [x2_p(1)/x2_p(3), x2_p(2)/x2_p(3)];
        % If the transformed point is inside the source image, save the
        % point coordinates
        if (x_e(1) > 1 && x_e(2) > 1 && x_e(1) < n && x_e(2) < m)
            indices(i,j) = 1;
            x2_aux_1(i,j) = x_e(1);
            x2_aux_2(i,j) = x_e(2);
=======
        x2_p = H\[x2_e;1]; %A\b for inv(a)*b
        x_e = [x2_p(1)/x2_p(3), x2_p(2)/x2_p(3)]; 
        if (x_e(1) > 1 && x_e(2) > 1 && x_e(1) < n && x_e(2) < m)
            Xq(i,j) = x_e(1);
            Yq(i,j) = x_e(2);
>>>>>>> origin/master
        end
    end
end

<<<<<<< HEAD
% Initialize interpolated image
I_interp = zeros(numel(I2(:,:,1)),3);

% Coordinates of the source image
[X,Y]=meshgrid(1:m,1:n);

% Interpolate each channel
I_interp(indices>0,1) = interp2(X,Y, double(I(:,:,1)),x2_aux_2(indices>0),x2_aux_1(indices>0));
I_interp(indices>0,2) = interp2(X,Y, double(I(:,:,2)),x2_aux_2(indices>0),x2_aux_1(indices>0));
I_interp(indices>0,3) = interp2(X,Y, double(I(:,:,3)),x2_aux_2(indices>0),x2_aux_1(indices>0));

% Reshape the interpolated image
I2 = reshape(I_interp,n2,m2,3);
=======
for k = 1:c
    I2(:,:,k) = interp2(single(I(:,:,k)), Xq, Yq);
end
I2(isnan(I2)) = 0;
>>>>>>> origin/master

end

