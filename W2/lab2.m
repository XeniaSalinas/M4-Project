%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 2: Image mosaics

addpath('sift');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute image correspondences

%% Open images

% imargb = imread('Data/llanes/llanes_a.jpg');
% imbrgb = imread('Data/llanes/llanes_b.jpg');
% imcrgb = imread('Data/llanes/llanes_c.jpg');

imargb = imread('Data/castle_int/0016_s.png');
imbrgb = imread('Data/castle_int/0015_s.png');
imcrgb = imread('Data/castle_int/0014_s.png');

% imargb = imread('Data/aerial/site13/frame00000.png');
% imbrgb = imread('Data/aerial/site13/frame00002.png');
% imcrgb = imread('Data/aerial/site13/frame00003.png');

ima = sum(double(imargb), 3) / 3 / 255;
imb = sum(double(imbrgb), 3) / 3 / 255;
imc = sum(double(imcrgb), 3) / 3 / 255;
sift_threshold = 0.01;

% imargb = imread('Data/aerial/site22/frame_00001.tif');
% imbrgb = imread('Data/aerial/site22/frame_00018.tif');
% imcrgb = imread('Data/aerial/site22/frame_00030.tif');
% ima = double(imargb) / 255.;
% imb = double(imbrgb) / 255.;
% imc = double(imcrgb) / 255.;
% sift_threshold = 0.02;

%% Compute SIFT keypoints
[points_a, desc_a] = sift(ima, 'Threshold', sift_threshold);
[points_b, desc_b] = sift(imb, 'Threshold', sift_threshold);
[points_c, desc_c] = sift(imc, 'Threshold', sift_threshold);

figure;
imshow(imargb);
hold on;
plot(points_a(1,:), points_a(2,:),'+y');
title('Image A keypoints');
hold off;

figure;
imshow(imbrgb); %image(imbrgb);
hold on;
plot(points_b(1,:), points_b(2,:),'+y');
title('Image B keypoints');
hold off;

figure;
imshow(imcrgb); %image(imcrgb);
hold on;
plot(points_c(1,:), points_c(2,:),'+y');
title('Image C keypoints');
hold off;

%% Match SIFT keypoints 

% between a and b
matches_ab = siftmatch(desc_a, desc_b);
figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), matches_ab, 'Stacking', 'v');

% between b and c
matches_bc = siftmatch(desc_b, desc_c);
figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), matches_bc, 'Stacking', 'v');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Compute the homography (DLT algorithm) between image pairs

%% Compute homography (normalized DLT) between a and b, play with the homography
th = 3;
xab_a = transpose(cart2hom(transpose(points_a(1:2, matches_ab(1,:)))));
xab_b = transpose(cart2hom(transpose(points_b(1:2, matches_ab(2,:)))));
[Hab, inliers_ab] = ransac_homography_adaptive_loop(xab_a, xab_b, th, 1000);

figure;
plotmatches(ima, imb, points_a(1:2,:), points_b(1:2,:), ...
    matches_ab(:,inliers_ab), 'Stacking', 'v');

vgg_gui_H(imargb, imbrgb, Hab);


%% Compute homography (normalized DLT) between b and c, play with the homography
xbc_b = transpose(cart2hom(transpose(points_b(1:2, matches_bc(1,:)))));
xbc_c = transpose(cart2hom(transpose(points_c(1:2, matches_bc(2,:)))));
[Hbc, inliers_bc] = ransac_homography_adaptive_loop(xbc_b, xbc_c, th, 1000); 

figure;
plotmatches(imb, imc, points_b(1:2,:), points_c(1:2,:), ...
    matches_bc(:,inliers_bc), 'Stacking', 'v');

vgg_gui_H(imbrgb, imcrgb, Hbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Build the mosaic

corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners);      % leave B unchanged   
iwa = apply_H_v2(imargb, Hab, corners);         % apply homography A --> B
iwc = apply_H_v2(imcrgb, inv(Hbc), corners);    % apply homography C --> B

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

% ToDo: compute the mosaic with castle_int images
% ToDo: compute the mosaic with aerial images set 13
% ToDo: compute the mosaic with aerial images set 22
% ToDo: comment the results in every of the four cases: say why it works or
%       does not work

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Refine the homography with the Gold Standard algorithm

%% Homography ab

% Set the non-homogeneous point coordinates of the point correspondences 
% we will refine with the geometric method
x = xab_a(1:2,inliers_ab);  
xp = xab_b(1:2,inliers_ab); 

Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hab(:) ; x(:) ];      % The parameters or independent variables

% Create a function that we pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares 
% E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
Y_initial = gs_errfunction( P0, Xobs ); 
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hab_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% Show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);


%% See differences in the keypoint locations

% Compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm
% xhat = P * x;
% xhatp = P * xp;
xhat = P(10:end);
len_x_hat = size(xhat,1)/2;
xhat = reshape(xhat, [2,len_x_hat]);
xhat = [xhat; ones(1,len_x_hat)];
xhatp = Hab_r * xhat;
xhat = euclid(xhat);
xhatp = euclid(xhatp);

figure;
imshow(imargb);%image(imargb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');
title('Image A, original (yellow) vs refined (blue) correspondences'); 
hold off;

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');
title('Image B, original (yellow) vs refined (blue) correspondences'); 
hold off;

%%  Homography bc

% Set the non-homogeneous point coordinates of the point correspondences 
% we will refine with the geometric method
x = xbc_b(1:2,inliers_bc);  
xp = xbc_c(1:2,inliers_bc); 

Xobs = [ x(:) ; xp(:) ];     % The column vector of observed values (x and x')
P0 = [ Hbc(:) ; x(:) ];      % The parameters or independent variables

% Create a function that we pass to the lsqnonlin function
% NOTE: gs_errfunction should return E(X) and not the sum-of-squares 
% E=sum(E(X).^2)) that we want to minimize. 
% (E(X) is summed and squared implicitly in the lsqnonlin algorithm.) 
Y_initial = gs_errfunction( P0, Xobs );
err_initial = sum( sum( Y_initial.^2 ));

options = optimset('Algorithm', 'levenberg-marquardt');
P = lsqnonlin(@(t) gs_errfunction(t, Xobs), P0, [], [], options);

Hbc_r = reshape( P(1:9), 3, 3 );
f = gs_errfunction( P, Xobs ); % lsqnonlin does not return f
err_final = sum( sum( f.^2 ));

% Show the geometric error before and after the refinement
fprintf(1, 'Gold standard reproj error initial %f, final %f\n', err_initial, err_final);


%% See differences in the keypoint locations

% ToDo: compute the points xhat and xhatp which are the correspondences
% returned by the refinement with the Gold Standard algorithm
% xhat = P * x;
% xhatp = P * xp;
xhat = P(10:end);
len_x_hat = size(xhat,1)/2;
xhat = reshape(xhat, [2,len_x_hat]);
xhat = [xhat; ones(1,len_x_hat)];
xhatp = Hbc_r * xhat;
xhat = euclid(xhat);
xhatp = euclid(xhatp);

figure;
imshow(imbrgb);%image(imbrgb);
hold on;
plot(x(1,:), x(2,:),'+y');
plot(xhat(1,:), xhat(2,:),'+c');
title('Image B, original (yellow) vs refined (blue) correspondences'); 
hold off;

figure;
imshow(imcrgb);%image(imcrgb);
hold on;
plot(xp(1,:), xp(2,:),'+y');
plot(xhatp(1,:), xhatp(2,:),'+c');
title('Image C, original (yellow) vs refined (blue) correspondences'); 
hold off;

%% Build mosaic
corners = [-400 1200 -100 650];
iwb = apply_H_v2(imbrgb, eye(3), corners); 
iwa = apply_H_v2(imargb, Hab_r, corners); 
iwc = apply_H_v2(imcrgb, inv(Hbc_r), corners); 

figure;
imshow(max(iwc, max(iwb, iwa)));%image(max(iwc, max(iwb, iwa)));axis off;
title('Mosaic A-B-C');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Calibration with a planar pattern

clear all;
%% Read template and images.
T     = imread('Data/calib/template.jpg');
I{1}  = imread('Data/calib/graffiti1.tif');
I{2}  = imread('Data/calib/graffiti2.tif');
I{3}  = imread('Data/calib/graffiti3.tif');
%I{4}  = imread('Data/calib/graffiti4.tif');
%I{5}  = imread('Data/calib/graffiti5.tif');
Tg = sum(double(T), 3) / 3 / 255;
Ig{1} = sum(double(I{1}), 3) / 3 / 255;
Ig{2} = sum(double(I{2}), 3) / 3 / 255;
Ig{3} = sum(double(I{3}), 3) / 3 / 255;

N = length(I);

%% Compute keypoints.
fprintf('Computing sift points in template... ');
[pointsT, descrT] = sift(Tg, 'Threshold', 0.05);
fprintf(' done\n');

points = cell(N,1);
descr = cell(N,1);
for i = 1:N
    fprintf('Computing sift points in image %d... ', i);
    [points{i}, descr{i}] = sift(Ig{i}, 'Threshold', 0.05);
    fprintf(' done\n');
end

%% Match and compute homographies.
H = cell(N,1);
for i = 1:N
    % Match against template descriptors.
    fprintf('Matching image %d... ', i);
    matches = siftmatch(descrT, descr{i});
    fprintf('done\n');

    % Fit homography and remove outliers.
    x1 = pointsT(1:2, matches(1, :));
    x2 = points{i}(1:2, matches(2, :));
    H{i} = 0;
        
    x1_h = transpose(cart2hom(transpose(x1)));
    x2_h = transpose(cart2hom(transpose(x2)));
    [H{i}, inliers] =  ransac_homography_adaptive_loop(x1_h, x2_h, 3, 1000);

    % Plot inliers.
    figure;
    plotmatches(Tg, Ig{i}, pointsT(1:2,:), points{i}(1:2,:), matches(:, inliers));

    % Play with the homography
%     vgg_gui_H(T, I{i}, H{i});
end

%% Compute the Image of the Absolute Conic
V = zeros(2*N,6); 

for n = 1:N
    % Compute row vector V11 (i=1, j=1)
    v11_t = compute_row_absolute_conic(H{n}, 1, 1);
    % Compute row vector V12 (i=1, j=2)
    v12_t = compute_row_absolute_conic(H{n}, 1, 2);
    % Compute row vector V22 (i=2, j=2)
    v22_t = compute_row_absolute_conic(H{n}, 2, 2);
    
    % 2 equations for homography n
    row = 2*(n-1)+1;
    V(row,:) = v12_t;
    V(row+1,:) = v11_t - v22_t;
end

[U,D,Uhat] = svd(V); 
% the solution is the last column of U_hat given by the SVD decomposition of V
omega = Uhat(:,end);  
w =[omega(1) omega(2) omega(3);
    omega(2) omega(4) omega(5);
    omega(3) omega(5) omega(6)];
 
%% Recover the camera calibration.

K = chol(w); 
    
% ToDo: in the report make some comments related to the obtained internal
%       camera parameters and also comment their relation to the image size

%% Compute camera position and orientation.
R = cell(N,1);
t = cell(N,1);
P = cell(N,1);
figure;hold;
for i = 1:N
    % ToDo: compute r1, r2, and t{i}
    r1 = K\H{1}(:,1);
    r2 = K\H{1}(:,2);
    t{i} = K\H{1}(:,3);
    
    % Solve the scale ambiguity by forcing r1 and r2 to be unit vectors.
    s = sqrt(norm(r1) * norm(r2)) * sign(t{i}(3));
    r1 = r1 / s;
    r2 = r2 / s;
    t{i} = t{i} / s;
    R{i} = [r1, r2, cross(r1,r2)];
    
    % Ensure R is a rotation matrix
    [U S V] = svd(R{i});
    R{i} = U * eye(3) * V';
   
    P{i} = K * [R{i} t{i}];
    plot_camera(P{i}, 800, 600, 200);
end

% ToDo: in the report explain how the optical center is computed in the
%       provided code

[ny,nx] = size(T);
p1 = [0 0 0]';
p2 = [nx 0 0]';
p3 = [nx ny 0]';
p4 = [0 ny 0]';
% Draw planar pattern
vgg_scatter_plot([p1 p2 p3 p4 p1], 'g');
% Paint image texture
surface('XData',[0 nx; 0 nx],'YData',[0 0; 0 0],'ZData',[0 0; -ny -ny],'CData',T,'FaceColor','texturemap');
colormap(gray);
axis equal;

%% Plot a static camera with moving calibration pattern.
figure; hold;
plot_camera(K * eye(3,4), 800, 600, 200);
% ToDo: complete the call to the following function with the proper
%       coordinates of the image corners in the new reference system
for i = 1:N
    vgg_scatter_plot( [K*p1   K*p2   K*p3   K*p4   K*p1], 'r');
end

%% Augmented reality: Plot some 3D points on every camera.
[Th, Tw] = size(Tg);
cube = [0 0 0; 1 0 0; 1 0 0; 1 1 0; 1 1 0; 0 1 0; 0 1 0; 0 0 0; 0 0 1; 1 0 1; 1 0 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 0 1; 0 0 0; 1 0 0; 1 0 0; 1 0 1; 1 0 1; 0 0 1; 0 0 1; 0 0 0; 0 1 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 0; 0 0 0; 0 1 0; 0 1 0; 0 1 1; 0 1 1; 0 0 1; 0 0 1; 0 0 0; 1 0 0; 1 1 0; 1 1 0; 1 1 1; 1 1 1; 1 0 1; 1 0 1; 1 0 0 ]';

X = (cube - .5) * Tw / 4 + repmat([Tw / 2; Th / 2; -Tw / 8], 1, length(cube));

for i = 1:N
    figure; colormap(gray);
    imagesc(Ig{i});
    hold on;
    x = euclid(P{i} * homog(X));
    vgg_scatter_plot(x, 'g');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. OPTIONAL: Detect the UPF logo in the two UPF images using the 
%%              DLT algorithm (folder "logos").
%%              Interpret and comment the results.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. OPTIONAL: Replace the logo of the UPF by the master logo
%%              in one of the previous images using the DLT algorithm.



