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
