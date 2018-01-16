%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 3: The geometry of two views 
% (application: photo-sequencing)

addpath('sift'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Compute the fundamental matrix

% Two camera matrices for testing purposes
P1 = eye(3,4);
c = cosd(15); s = sind(15);
R = [c -s 0; s c 0; 0 0 1];
t = [.3 0.1 0.2]';
P2 = [R t];
n = 8;
X = [rand(3,n); ones(1,n)] + [zeros(2,n); 3 * ones(1,n); zeros(1,n)];
x1_test = P1 * X;
x2_test = P2 * X;

% Estimated fundamental matrix
% ToDo: create the following function that estimates F using the normalised 8 point algorithm
F_es = fundamental_matrix(x1_test, x2_test);

% Real fundamental matrix
% K and K_ should be I as P1 = K[I | 0] = I and P2 = K_[R | t] = [R | t], 
% but let's compute them:
% get K and K' from the QR decomposition of P1 and P2:
[~, K] = qr(P1(:,1:end-1));
[~, K_] = qr(P2(:,1:end-1));
%enforce that K and K_ have positive diagonal entries
D = diag(sign(diag(K)));
D_ = diag(sign(diag(K_)));
K = (D*K);
K_ = (D_*K_);
%Compute the matrix Tx:
Tx = [0, -t(3), t(2); t(3), 0, -t(1);-t(2), t(1), 0];

F_gt = inv(K_)'*(Tx*R)*inv(K); % ToDo: write the expression of the real fundamental matrix for P1 and P2

% Evaluation: these two (normalized) matrices should be very similar
norm_F_gt = F_gt / norm(F_gt);
norm_F_es = F_es / norm(F_es);

if norm(norm_F_gt - norm_F_es) < 0.1
    disp('Ground-truth and estimated fundamental matrix are similar');
else
    warning('Ground-truth and estimated fundamental matrix are different');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Robustly fit fundamental matrix

% Read images
im1rgb = imread('Data/0000_s.png');
im2rgb = imread('Data/0001_s.png');
im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;

% show images
figure(1);
subplot(1,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(1,2,2); imshow(im2rgb); axis image; title('Image 2');


%% Compute SIFT keypoints

% (make sure that the sift folder provided in lab2 is on the path)

[points_1, desc_1] = sift(im1, 'Threshold', 0.01);
[points_2, desc_2] = sift(im2, 'Threshold', 0.01);

%% Match SIFT keypoints between a and b
matches = siftmatch(desc_1, desc_2);
figure(2);
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches, 'Stacking', 'v');

% p1 and p2 contain the homogeneous coordinates of the matches
p1 = [points_1(1:2, matches(1,:)); ones(1, length(matches))];
p2 = [points_2(1:2, matches(2,:)); ones(1, length(matches))];

% ToDo: create this function (you can use as a basis 'ransac_homography_adaptive_loop.m')
[F, inliers] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 

% show inliers
figure(3);
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches(:,inliers), 'Stacking', 'v');
title('Inliers');

vgg_gui_F(im1rgb, im2rgb, F');


%% Plot some epipolar lines

l2 = F*p1; % epipolar lines in image 2 
l1 = F'*p2; % epipolar lines in image 1

% choose three random indices
permuted_inliers = inliers(randperm(length(inliers)));
m1 = permuted_inliers(1);
m2 = permuted_inliers(2);
m3 = permuted_inliers(3);

% image 1 (plot the three points and their corresponding epipolar lines)
figure(4);
imshow(im1rgb);
hold on;
plot(p1(1, m1), p1(2, m1), '+g');
plot_homog_line(l1(:, m1), 'y');

plot(p1(1, m2), p1(2, m2), '+g');
plot_homog_line(l1(:, m2), 'r');

plot(p1(1, m3), p1(2, m3), '+g');
plot_homog_line(l1(:, m3), 'b');

title('Epipolar lines (image 1)');

% image 2 (plot the three points and their corresponding epipolar lines)
figure(5);
imshow(im2rgb);
hold on;
plot(p2(1, m1), p2(2, m1), '+g');
plot_homog_line(l2(:, m1), 'y');

plot(p2(1, m2), p2(2, m2), '+g');
plot_homog_line(l2(:, m2), 'r');

plot(p2(1, m3), p2(2, m3), '+g');
plot_homog_line(l2(:, m3), 'b');

title('Epipolar lines (image 2)');

%% Compare RANSAC with non-robust version of 8-point algorithm
F_no_robust = fundamental_matrix(p1, p2);

l2_no_robust = F_no_robust*p1;  % epipolar lines in image 2 using 8-point, no RANSAC
l1_no_robust = F_no_robust'*p2; % epipolar lines in image 1 using 8-point, no RANSAC

% image 1
figure(6);
imshow(im1rgb);
hold on;
plot(p1(1, m1), p1(2, m1), '+g');
plot_homog_line(l1_no_robust(:, m1), 'y');

plot(p1(1, m2), p1(2, m2), '+g');
plot_homog_line(l1_no_robust(:, m2), 'r');

plot(p1(1, m3), p1(2, m3), '+g');
plot_homog_line(l1_no_robust(:, m3), 'b');

title('Epipolar lines, not RANSAC (image 1)');

% image 2
figure(7);
imshow(im2rgb);
hold on;
plot(p2(1, m1), p2(2, m1), '+g');
plot_homog_line(l2_no_robust(:, m1), 'y');

plot(p2(1, m2), p2(2, m2), '+g');
plot_homog_line(l2_no_robust(:, m2), 'r');

plot(p2(1, m3), p2(2, m3), '+g');
plot_homog_line(l2_no_robust(:, m3), 'b');

title('Epipolar lines, not RANSAC (image 2)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Photo-sequencing with aerial images

% In this part we will compute a simplified version of the algorithm
% explained in the Photo-sequencing paper. 
% Since we do not have two images
% taken from roughly the same viewpoint at two different time instants we
% will manually pick a dynamic point corresponding to a point in a van 
% (identified by index 'idx_car_I1') and the projection of its 3D trajectory 
% in the reference image. Then we will compute the projection (to the reference image) 
% of three points on this 3D trajectory at three different time instants 
% (corresponding to the time when the three other provided images where taken). 

clear all;

% Read images
im1rgb = imread('Data/frame_00000.tif');
im2rgb = imread('Data/frame_00001.tif');
im3rgb = imread('Data/frame_00002.tif');
im4rgb = imread('Data/frame_00003.tif');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure(8);
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

%% Compute the fundamental matrices

% Take image im1 as reference image (image 1) and compute the fundamental 
% matrices needed for computing the trajectory of point idx_car_I1
% (use the SIFT keypoints previously computed)
matches_2 = siftmatch(desc_1, desc_2);
p1 = [points_1(1:2, matches_2(1,:)); ones(1, length(matches_2))];
p2 = [points_2(1:2, matches_2(2,:)); ones(1, length(matches_2))];
[F_2, inliers_2] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 
% show inliers
figure(9);
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches_2(:,inliers_2), 'Stacking', 'v');
title('Inliers 1 and 2');

matches_3 = siftmatch(desc_1, desc_3);
p1 = [points_1(1:2, matches_3(1,:)); ones(1, length(matches_3))];
p3 = [points_3(1:2, matches_3(2,:)); ones(1, length(matches_3))];
[F_3, inliers_3] = ransac_fundamental_matrix(p1, p3, 2.0, 1000);
% show inliers
figure(10);
plotmatches(im1, im3, points_1(1:2,:), points_3(1:2,:), matches_3(:,inliers_3), 'Stacking', 'v');
title('Inliers 1 and 3');

matches_4 = siftmatch(desc_1, desc_4);
p1 = [points_1(1:2, matches_4(1,:)); ones(1, length(matches_4))];
p4 = [points_4(1:2, matches_4(2,:)); ones(1, length(matches_4))];
[F_4, inliers_4] = ransac_fundamental_matrix(p1, p4, 2.0, 1000); 
% show inliers
figure(11);
plotmatches(im1, im4, points_1(1:2,:), points_4(1:2,:), matches_4(:,inliers_4), 'Stacking', 'v');
title('Inliers 1 and 4');


%% Plot the car trajectory (keypoint idx_car_I1 in image 1)

idx_car_I1 = 1197;
idx_car_I2 = matches_2(2,matches_2(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 2
idx_car_I3 = matches_3(2,matches_3(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 3
idx_car_I4 = matches_4(2,matches_4(1,:)==idx_car_I1); % ToDo: identify the corresponding point of idx_car_I1 in image 4

% coordinates (in image 1) of the keypoint idx_car_I1 (point in a van). 
% point1_1 is the projection of a 3D point in the 3D trajectory of the van
point1_1 = [points_1(1:2,idx_car_I1)' 1]';
% coordinates (in image 1) of another 3D point in the same 3D trajectory of
% the van
point1_2 = [334 697 1]'; % (this is a given data)

% l1 is the projection of the 3D trajectory of keypoint idx_car_I1
% (it is the line that joins point1_1 and point1_2)
l1 = cross(point1_1, point1_2);% ToDo: compute the line
% plot the line
figure(12);
imshow(im1);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y', 'DisplayName', 'Van trajectory');
plot(points_1(1, idx_car_I1), points_1(2, idx_car_I1), 'y*', 'DisplayName', 'Van (image 1)');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I2 in image 2
point2 = [points_2(1:2,idx_car_I2)' 1]';
% ToDo: compute the epipolar line of point2 in the reference image
l2 = F_2' * point2;
% plot the epipolar line
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'c', 'DisplayName', 'Epipolar line (image 2)');
% ToDo: compute the projection of point idx_car_I2 in the reference image 
pi2 = cross(l1, l2);
% plot this point
plot(pi2(1)/pi2(3), pi2(2)/pi2(3), 'c*', 'DisplayName', 'Van (image 2)');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I3 in image 3
point3 = [points_3(1:2,idx_car_I3)' 1]';
% ToDo: compute the epipolar line of point3 in the reference image
l3 = F_3' * point3;
% plot the epipolar line
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b', 'DisplayName', 'Epipolar line (image 3)');
% ToDo: compute the projection of point idx_car_I3 in the reference image
pi3 = cross(l1, l3);
plot(pi3(1)/pi3(3), pi3(2)/pi3(3), 'b*', 'DisplayName', 'Van (image 3)');

% ToDo: write the homogeneous coordinates of the corresponding point of idx_car_I4 in image 4
point4 = [points_4(1:2,idx_car_I4)' 1]';
% ToDo: compute the epipolar line of point4 in the reference image
l4 = F_4' * point4;
% plot the epipolar line
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g', 'DisplayName', 'Epipolar line (image 4)');
% ToDo: compute the projection of point idx_car_I4 in the reference image
pi4 = cross(l1, l4);
plot(pi4(1)/pi4(3), pi4(2)/pi4(3), 'g*', 'DisplayName', 'Van (image 4)');

% Plot legend
legend('show');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. OPTIONAL: Photo-sequencing with your own images

% 4.1 Take a set of images of a moving scene from different viewpoints at 
%     different time instants. At least two images have to be taken from
%     roughly the same location by the same camera.
clear all;

% Read images
im1rgb = imread('Data/Photo_seq/IMG_20180116_122901.jpg'); %Reference image
im2rgb = imread('Data/Photo_seq/IMG_20180116_122909.jpg'); %Same point of view
im3rgb = imread('Data/Photo_seq/IMG_20180116_122927.jpg');
im4rgb = imread('Data/Photo_seq/IMG_20180116_122944.jpg');

im1 = sum(double(im1rgb), 3) / 3 / 255;
im2 = sum(double(im2rgb), 3) / 3 / 255;
im3 = sum(double(im3rgb), 3) / 3 / 255;
im4 = sum(double(im4rgb), 3) / 3 / 255;

% show images
figure(13);
subplot(2,2,1); imshow(im1rgb); axis image; title('Image 1');
subplot(2,2,2); imshow(im2rgb); axis image; title('Image 2');
subplot(2,2,3); imshow(im3rgb); axis image; title('Image 3');
subplot(2,2,4); imshow(im4rgb); axis image; title('Image 4');

% 4.2 Implement the first part (until line 16) of the Algorithm 1 of the 
%     Photo-sequencing paper with a selection of the detected dynamic
%     features. You may reuse the code generated for the previous question.
% Compute SIFT keypoints
[points_1, desc_1] = sift(im1, 'Threshold', 0.015); % Do not change this threshold!
[points_2, desc_2] = sift(im2, 'Threshold', 0.015);
[points_3, desc_3] = sift(im3, 'Threshold', 0.015);
[points_4, desc_4] = sift(im4, 'Threshold', 0.015);

%% Compute the fundamental matrices

% Take image im1 as reference image (image 1) and compute the fundamental 
% matrices needed for computing the trajectory of point idx_car_I1
% (use the SIFT keypoints previously computed)
matches_2 = siftmatch(desc_1, desc_2);
p1 = [points_1(1:2, matches_2(1,:)); ones(1, length(matches_2))];
p2 = [points_2(1:2, matches_2(2,:)); ones(1, length(matches_2))];
[F_2, inliers_2] = ransac_fundamental_matrix(p1, p2, 2.0, 1000); 
% show inliers
figure(14);
plotmatches(im1, im2, points_1(1:2,:), points_2(1:2,:), matches_2(:,:), 'Stacking', 'v');
title('Inliers 1 and 2');

matches_3 = siftmatch(desc_1, desc_3);
p1 = [points_1(1:2, matches_3(1,:)); ones(1, length(matches_3))];
p3 = [points_3(1:2, matches_3(2,:)); ones(1, length(matches_3))];
[F_3, inliers_3] = ransac_fundamental_matrix(p1, p3, 2.0, 1000);
% show inliers
figure(15);
plotmatches(im1, im3, points_1(1:2,:), points_3(1:2,:), matches_3(:,:), 'Stacking', 'v');
title('Inliers 1 and 3');

matches_4 = siftmatch(desc_1, desc_4);
p1 = [points_1(1:2, matches_4(1,:)); ones(1, length(matches_4))];
p4 = [points_4(1:2, matches_4(2,:)); ones(1, length(matches_4))];
[F_4, inliers_4] = ransac_fundamental_matrix(p1, p4, 2.0, 1000); 
% show inliers
figure(16);
plotmatches(im1, im4, points_1(1:2,:), points_4(1:2,:), matches_4(:,:), 'Stacking', 'v');
title('Inliers 1 and 4');

%% Compute the dynamic points
% Error between estimates and observations
dist = sqrt((points_1(1,matches_2(1,:))-points_2(1,matches_2(2,:))).^2 + ...
    (points_1(2,matches_2(1,:))-points_2(2,matches_2(2,:))).^2);
dynamic = dist > mean(dist);
%show points
figure(17);
imshow(im1);
hold on;
plot(points_1(1,matches_2(1,:)),points_1(2,matches_2(1,:)), 'y*','DisplayName', 'Points in image 1');
plot(points_2(1,matches_2(2,not(dynamic))), points_2(2,matches_2(2,not(dynamic))), 'b*','DisplayName', 'Static points in image 2');
plot(points_2(1,matches_2(2,dynamic)), points_2(2,matches_2(2,dynamic)), 'r*','DisplayName', 'Dynamic points in image 2');
% Plot legend
legend('show');

%% Compute images order
for dyn = matches_2(1,dynamic)

    idx_doll_I1 = dyn;
    idx_doll_I2 = matches_2(2,matches_2(1,:)==idx_doll_I1); % ToDo: identify the corresponding point of idx_doll_I1 in image 2
    idx_doll_I3 = matches_3(2,matches_3(1,:)==idx_doll_I1); % ToDo: identify the corresponding point of idx_doll_I1 in image 3
    idx_doll_I4 = matches_4(2,matches_4(1,:)==idx_doll_I1); % ToDo: identify the corresponding point of idx_doll_I1 in image 4

    %Check if the dynamic point is in all the images
    if (~isempty(idx_doll_I2) && ~isempty(idx_doll_I3) && ~isempty(idx_doll_I4))
        % coordinates (in image 1) of the keypoint idx_car_I1 (point in a van). 
        % point1_1 is the projection of a 3D point in the 3D trajectory of the van
        point1_1 = [points_1(1:2,idx_doll_I1)' 1]';
        % coordinates (in image 1) of another 3D point in the same 3D trajectory of the doll
        point1_2 = [points_2(1:2,idx_doll_I2)' 1]';

        % l1 is the projection of the 3D trajectory of keypoint idx_doll_I1
        % (it is the line that joins point1_1 and point1_2)
        l1 = cross(point1_1, point1_2);% ToDo: compute the line

        % Compute the projection of point idx_doll_I1 in the reference image 
        pi1 = point1_1;
        alpha1 = (pi1(1) - point1_1(1))./ (point1_1(1) - point1_2(1));

        % Write the homogeneous coordinates of the corresponding point of idx_doll_I2 in image 2
        point2 = [points_2(1:2,idx_doll_I2)' 1]';
        % Compute the epipolar line of point2 in the reference image
        l2 = F_2' * point2;
        % Compute the projection of point idx_doll_I2 in the reference image 
        pi2 = cross(l1, l2);
        alpha2 = (pi2(1) - point1_1(1))./ (point1_1(1) - point1_2(1));

        % Write the homogeneous coordinates of the corresponding point of idx_doll_I3 in image 3
        point3 = [points_3(1:2,idx_doll_I3)' 1]';
        % Compute the epipolar line of point3 in the reference image
        l3 = F_3' * point3;
        % Compute the projection of point idx_doll_I3 in the reference image 
        pi3 = cross(l1, l3);
        alpha3 = (pi3(1) - point1_1(1))./ (point1_1(1) - point1_2(1));

        % Write the homogeneous coordinates of the corresponding point of idx_doll_I4 in image 4
        point4 = [points_4(1:2,idx_doll_I4)' 1]';
        % Compute the epipolar line of point4 in the reference image
        l4 = F_4' * point4;
        % Compute the projection of point idx_doll_I4 in the reference image 
        pi4 = cross(l1, l4);
        alpha4 = (pi4(1) - point1_1(1))./ (point1_1(1) - point1_2(1));

        [sigma, order] = sort([alpha1 alpha2, alpha3, alpha4]);
        order
    end

end
