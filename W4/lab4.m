%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 4: Reconstruction from two views (knowing internal camera parameters) 


addpath('sift'); % ToDo: change 'sift' to the correct path where you have the sift functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Triangulation

% ToDo: create the function triangulate.m that performs a triangulation
%       with the homogeneous algebraic method (DLT)
%
%       The entries are (x1, x2, P1, P2, imsize), where:
%           - x1, and x2 are the Euclidean coordinates of two matching 
%             points in two different images.
%           - P1 and P2 are the two camera matrices
%           - imsize is a two-dimensional vector with the image size

%% Test the triangulate function
% Use this code to validate that the function triangulate works properly

P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X_test = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = euclid(P1 * X_test);
x2_test = euclid(P2 * X_test);

N_test = size(x1_test,2);
X_trian = zeros(4,N_test);
for i = 1:N_test
    X_trian(:,i) = triangulate(x1_test(:,i), x2_test(:,i), P1, P2, [2 2]);
end

% error
triangulation_error = euclid(X_test) - euclid(X_trian);
display(triangulation_error);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Reconstruction from two views

%% Read images
Irgb{1} = imread('Data/0001_s.png');
Irgb{2} = imread('Data/0002_s.png');
I{1} = sum(double(Irgb{1}), 3) / 3 / 255;
I{2} = sum(double(Irgb{2}), 3) / 3 / 255;
[h,w] = size(I{1});


%% Compute keypoints and matches.
points = cell(2,1);
descr = cell(2,1);
for i = 1:2
    [points{i}, descr{i}] = sift(I{i}, 'Threshold', 0.01);
    points{i} = points{i}(1:2,:);
end

matches = siftmatch(descr{1}, descr{2});

% Plot matches.
figure();
plotmatches(I{1}, I{2}, points{1}, points{2}, matches, 'Stacking', 'v');


%% Fit Fundamental matrix and remove outliers.
x1 = points{1}(:, matches(1, :));
x2 = points{2}(:, matches(2, :));
[F, inliers] = ransac_fundamental_matrix(homog(x1), homog(x2), 2.0);

% Plot inliers.
inlier_matches = matches(:, inliers);
figure;
plotmatches(I{1}, I{2}, points{1}, points{2}, inlier_matches, 'Stacking', 'v');

x1 = points{1}(:, inlier_matches(1, :));
x2 = points{2}(:, inlier_matches(2, :));

vgg_gui_F(Irgb{1}, Irgb{2}, F');


%% Compute candidate camera matrices.

% Camera calibration matrix
K = [2362.12 0 1520.69; 0 2366.12 1006.81; 0 0 1];
scale = 0.3;
H = [scale 0 0; 0 scale 0; 0 0 1];
K = H * K;

% ToDo: Compute the Essential matrix from the Fundamental matrix
E = transpose(K) * F * K;
P1 = K*eye(3,4);
%Compute the extrinsic parameters of P2 from the Essential matrix
[U,D,V]=svd(E);

W=[0,-1,0; 1,0,0;0,0,1];
Z=[0,1,0;-1,0,0;0,0,0];

% Direct way to compute the (unsigned) translation vector --> last vector
% of U
t=U(:,end);

% Indirect way using the SVD decomposition of the skew-symmetric matrix
% [Tx] --> SVD([Tx]) --> last column of V
Tx_skew_symmetric = U * (- Z) * transpose(U);
[~, ~, V_Tx] = svd(Tx_skew_symmetric);
t_prima = V(:, end);

% The translation vectors must be equal (possibly with changed signs)
display(t);
display(t_prima);

% Compute rotation
R1=U*W*V';
R2=U*W'*V';


% HINT: You may get improper rotations; in that case you need to change
%       their sign.
% Let R be a rotation matrix, you may check:
% if det(R) < 0
%     R = -R;
% end
 if det(R1) < 0
     R1 = -R1;
 end
 if det(R2) < 0
     R2 = -R2;
 end

Pc2 = {};
Pc2{1} = K * [R1, t];
Pc2{2} = K * [R1, -t];
Pc2{3} = K * [R2 , t];
Pc2{4} = K * [R2 , -t];

% plot the first camera and the four possible solutions for the second
figure;
plot_camera(P1,w,h, 'k');
plot_camera(Pc2{1},w,h,'y');
plot_camera(Pc2{2},w,h,'r');
plot_camera(Pc2{3},w,h,'g');
plot_camera(Pc2{4},w,h,'b');

% create legend
hold on;
handler = zeros(5, 1);
handler(1) = plot(NaN,NaN,'-k');
handler(2) = plot(NaN,NaN,'-y');
handler(3) = plot(NaN,NaN,'-r');
handler(4) = plot(NaN,NaN,'-g');
handler(5) = plot(NaN,NaN,'-b');
legend(handler, 'Camera 1', 'Camera 2 (R1 | t)', 'Camera 2 (R1 | -t)', 'Camera 2 (R2 | t)', 'Camera 2 (R2 | -t)');
hold off;

%% Reconstruct structure
% ToDo: Choose a second camera candidate by triangulating a match.

% Choose random match, use seed for reproducibility
seed = 123456;
rng(seed);
match_idx = randi(length(x1));

selected_camera = 0;
P2 = 0;
for i=1:4
    X=triangulate(x1(:,match_idx), x2(:,match_idx), P1, Pc2{i}, [w h]);
    
    % Projection into image planes
    x_p1 = P1 * X;
    x_p2 = Pc2{i} * X;
    
    if x_p1(3) > 0 && x_p2(3) > 0
        P2 = Pc2{i};
        selected_camera = i;
        break;
    end
end
fprintf('Selected camera: %d\n', selected_camera);

% Triangulate all matches.
N = size(x1,2);
X = zeros(4,N);
for i = 1:N
    X(:,i) = triangulate(x1(:,i), x2(:,i), P1, P2, [w h]);
end



%% Plot with colors
r = interp2(double(Irgb{1}(:,:,1)), x1(1,:), x1(2,:));
g = interp2(double(Irgb{1}(:,:,2)), x1(1,:), x1(2,:));
b = interp2(double(Irgb{1}(:,:,3)), x1(1,:), x1(2,:));
Xe = euclid(X);
figure; hold on;
% Plot real cameras
plot_camera(P1,w,h,'b');
plot_camera(P2,w,h, 'r');
% Plot alternative cameras as well, to see how they would be positioned
colors = ['k', 'y', 'g'];
color_idx = 1;
for i = 1:4
    if i == selected_camera
        continue
    end
    plot_camera(Pc2{i},w,h, colors(color_idx));
    color_idx = color_idx + 1;
end
% Plot 3D points
for i = 1:length(Xe)
    scatter3(Xe(1,i), Xe(3,i), -Xe(2,i), 5^2, [r(i) g(i) b(i)]/255, 'filled');
end
axis equal;
% Create legend
hold on;
handler = zeros(5, 1);
handler(1) = plot(NaN,NaN,'-b');
handler(2) = plot(NaN,NaN,'-r');
handler(3) = plot(NaN,NaN,'-k');
handler(4) = plot(NaN,NaN,'-y');
handler(5) = plot(NaN,NaN,'-g');
legend(handler, 'Camera 1', 'Camera 2', 'Alt. Camera 2', 'Alt. Camera 2', 'Alt. Camera 2');
hold off;


%% Compute reprojection error.

% ToDo: compute the reprojection errors
%       plot the histogram of reprojection errors, and
%       plot the mean reprojection error
x1h=homog(x1);
x2h=homog(x2);

x1p=P1*X;
x2p=P2*X;

%Normalize projected points
x2p=x2p./x2p(3,:);
x1p=x1p./x1p(3,:);

re=sum((x1h-x1p).^2) + sum((x2h-x2p).^2); % Reprojection error

%Histogram of reprojection errors
figure(); hist(re, 50);

%Mean reprojection error
re_mean=mean(re);
hold on
plot([re_mean,re_mean],ylim,'r')
hold off

disp(['Mean reprojection error: ' num2str(re_mean)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Depth map computation with local methods (SSD)

% Data images: 'scene1.row3.col3.ppm','scene1.row3.col4.ppm'
% Disparity ground truth: 'truedisp.row3.col3.pgm'

% Write a function called 'stereo_computation' that computes the disparity
% between a pair of rectified images using a local method based on a matching cost 
% between two local windows.
% 
% The input parameters are 5:
% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)
%
% In this part we ask to implement only the SSD cost
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.
%
% Note 1: Use grayscale images
% Note 2: For this first set of images use 0 as minimum disparity 
% and 16 as the the maximum one.

% Load images and gt disparity map
stereo_img_l = imread('Data/scene1.row3.col3.ppm');
stereo_img_r = imread('Data/scene1.row3.col4.ppm');
stereo_img_l = sum(double(stereo_img_l), 3) / 3 / 255;
stereo_img_r = sum(double(stereo_img_r), 3) / 3 / 255;

gt_disp_map = imread('Data/truedisp.row3.col3.pgm');

% Parameters
min_disparity = 0;
max_disparity = 16;
window_size = 3;
matching_cost = 'SSD';

% Compute disparity map
disp_map = stereo_computation( ...
    stereo_img_l, stereo_img_r, ...
    min_disparity, max_disparity, ...
    window_size, matching_cost ...
);

% Compare disparity maps
figure;
imshow(disp_map / max(max(disp_map)));
title('Computed disparity map');

figure;
imshow(gt_disp_map);
title('Ground-truth disparity map');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Depth map computation with local methods (NCC)

% Complete the previous function by adding the implementation of the NCC
% cost.
%
% Evaluate the results changing the window size (e.g. 3x3, 9x9, 20x20,
% 30x30) and the matching cost. Comment the results.

% Load images and gt disparity map
stereo_img_l = imread('Data/scene1.row3.col3.ppm');
stereo_img_r = imread('Data/scene1.row3.col4.ppm');
stereo_img_l = sum(double(stereo_img_l), 3) / 3 / 255;
stereo_img_r = sum(double(stereo_img_r), 3) / 3 / 255;

gt_disp_map = imread('Data/truedisp.row3.col3.pgm');

% Parameters
min_disparity = 0;
max_disparity = 16;
window_size = 30;
matching_cost = 'NCC';

% Compute disparity map
disp_map = stereo_computation( ...
    stereo_img_l, stereo_img_r, ...
    min_disparity, max_disparity, ...
    window_size, matching_cost ...
);

% Compare disparity maps
figure;
imshow(disp_map / max(max(disp_map)));
title('Computed disparity map');

figure;
imshow(gt_disp_map);
title('Ground-truth disparity map');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Depth map computation with local methods

% Data images: '0001_rectified_s.png','0002_rectified_s.png'

% Test the functions implemented in the previous section with the facade
% images. Try different matching costs and window sizes and comment the
% results.
% Notice that in this new data the minimum and maximum disparities may
% change.

% Load images and gt disparity map
stereo_img_l = imread('Data/0001_rectified_s.png');
stereo_img_r = imread('Data/0002_rectified_s.png');
stereo_img_l = sum(double(stereo_img_l), 3) / 3 / 255;
stereo_img_r = sum(double(stereo_img_r), 3) / 3 / 255;

% Parameters
min_disparity = 0;
max_disparity = 16;
window_size = 3;
matching_cost = 'SSD';

% Compute disparity map
disp_map = stereo_computation( ...
    stereo_img_l, stereo_img_r, ...
    min_disparity, max_disparity, ...
    window_size, matching_cost ...
);

% Compare disparity maps
figure;
imshow(disp_map / max(max(disp_map)));
title('Computed disparity map');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Bilateral weights

% Modify the 'stereo_computation' so that you can use bilateral weights (or
% adaptive support weights) in the matching cost of two windows.
% Reference paper: Yoon and Kweon, "Adaptive Support-Weight Approach for Correspondence Search", IEEE PAMI 2006
%
% Comment the results and compare them to the previous results (no weights).
%
% Note: Use grayscale images (the paper uses color images)

% Load images and gt disparity map
stereo_img_l = imread('Data/scene1.row3.col4.ppm');
stereo_img_r = imread('Data/scene1.row3.col3.ppm');
stereo_img_l = sum(double(stereo_img_l), 3) / 3 / 255;
stereo_img_r = sum(double(stereo_img_r), 3) / 3 / 255;

% Parameters
min_disparity = 0;
max_disparity = 16;
window_size = 20;
matching_cost = 'BILATERAL';

% Compute disparity map
disp_map = stereo_computation( ...
    stereo_img_l, stereo_img_r, ...
    min_disparity, max_disparity, ...
    window_size, matching_cost ...
);

% Compare disparity maps
figure;
imshow(disp_map / max(max(disp_map)));
title('Computed disparity map');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Stereo computation with Belief Propagation

% Use the UGM library used in module 2 and implement a  
% stereo computation method that minimizes a simple stereo energy with 
% belief propagation. 
% For example, use an L2 or L1 pixel-based data term (SSD or SAD) and 
% the same regularization term you used in module 2. 
% Or pick a stereo paper (based on belief propagation) from the literature 
% and implement it. Pick a simple method or just simplify the method they propose.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth computation with Plane Sweeping

% Implement the plane sweeping method explained in class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  Depth map fusion 

% In this task you are asked to implement the depth map fusion method
% presented in the following paper:
% B. Curless and M. Levoy. A Volumetric Method for Building Complex
% Models from Range Images. In Proc. SIGGRAPH, 1996.
%
% 1. Use the set of facade images 00xx_s.png to compute depth maps 
% corresponding to different views (and optionally from different pairs of 
% images for the same view).
% 2. Then convert each depth map to a signed distance function defined in 
% a disretized volume (using voxels).
% 3. Average the different signed distance functions, the resulting 
% signed distance is called D.
% 4. Set as occupied voxels (those representing the surface) those 
% where D is very close to zero. The rest of voxels will be considered as 
% empty.
%
% For that you need to compute a depth map from a pair of views in general
% position (non rectified). Thus, you may either use the plane sweep
% algorithm (if you did it) or the local method for estimating depth
% (mandatory task) together with the following rectification method which 
% has an online demo available: 
% http://demo.ipol.im/demo/m_quasi_euclidean_epipolar_rectification/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTIONAL:  New view synthesis

% In this task you are asked to implement part of the new view synthesis method
% presented in the following paper:
% S. Seitz, and C. Dyer, View morphing, Proc. ACM SIGGRAPH 1996.

% You will use a pair of rectified stereo images (no need for prewarping
% and postwarping stages) and their corresponding ground truth disparities
% (folder "new_view").
% Remember to take into account occlusions as explained in the lab session.
% Once done you can apply the code to the another pair of rectified images 
% provided in the material and use the estimated disparities with previous methods.
