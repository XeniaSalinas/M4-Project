close all, clear, clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lab 1: Image rectification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Applying image transformations

% Create the function  "apply_H" that gets as input a homography and
% an image and returns the image transformed by the homography.
% The size of the transformed image has to be automatically set so as to 
% contain the whole transformed image.
% At some point you will need to interpolate the image values at some points,
% you may use the Matlab function "interp2" for that.

I=imread('Data/0005_s.png'); % we have to be in the proper folder

%% 1.1. Similarities
% Generate a matrix H which produces a similarity transformation
s = 1;                    % scale factor
angle = 30;               % rotation angle
t = [-1, 4];              % translation 
H = H_similarity_matrix(s, angle, t);
% Get the transformed image by the homography
I2 = apply_H(I, H);

figure(1), imshow(I), title('Image 0005_s');
figure(2), imshow(uint8(I2)), title('Image 0005_s with similarity transformation');

%% 1.2. Affinities

% Generate a matrix H which produces an affine transformation
A = [1 -0.1;
     0.5 0.8];
T = [20; 30];
H = H_affine_matrix(A, T);
I2 = apply_H(I, H);
figure(3), imshow(uint8(I2));
title('Image 0005_s with affinity transformation');

%% 1.2.1 Affinity decomposition
% Decompose the affinity in four transformations: two
% rotations, a scale, and a translation
[U,D,V] = svd(A);
R_theta = U * transpose(V);
R_phi = transpose(V);
S = D(1:2,1:2);

% Verify that the product of the four previous transformations
% produces the same matrix H as above
tolerance = 1e-4;
A2 = R_theta * transpose(R_phi) * S * R_phi;
H_decomposed = H_affine_matrix(A2, T);
if abs(sum(sum(H - H_decomposed))) > tolerance
    error('H is no equal to its decomposition');
end

% Verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I3 = apply_H(I, [eye(2) T; 0 0 1]);
I3 = apply_H(I3,[R_phi T; 0 0 1]);
I3 = apply_H(I3,[S [0; 0]; 0 0 1]);
I3 = apply_H(I3,[transpose(R_phi) [0; 0]; 0 0 1]);
I3 = crop_transformed_image(I3);
I3 = apply_H(I3,[R_theta [0; 0]; 0 0 1]);
I3 = crop_transformed_image(I3);
figure(4); 
imshow(uint8(I3));
title('Image 0005_s with affinity (decomposed) transformation');
% Compare the image transformed using H and the image transformed using the
% decomposition of H.
[rows, cols, c] = size(I2);
I3_resize = imresize(I3, [rows, cols]);
error_im = uint8(I2 - I3_resize);
mae = mean2(error_im);  % mean absolute difference
fprintf('Mean average error: %.3f\n', mae);
figure(12); imshow(error_im); 

%% 1.3 Projective transformations (homographies)

% Generate a matrix H which produces a projective transformation
A = [0.9 0.12;
     0.05 1.1];
t = [3; -4];
v = [0.0005 0.001];
H_projective = H_projective_matrix(A, t, v);
I_projective = apply_H(I, H_projective);
figure(5); 
imshow(uint8(I_projective));
title('Image 0005_s with projective transformation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification

I=imread('Data/0000_s.png');
line_info_file = load('Data/0000_s_info_lines.txt');
% Choose the image points
[H_affine, I_affine] = affine_rectification(I,line_info_file,424,240,712,565);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

[H_metric, I_metric] = metric_rectification(I_affine,H_affine,line_info_file,424,240,712,565);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

I=imread('Data/0001_s.png');
I = I(1:614, 1:533, :);
line_info_file = load('Data/0001_s_info_lines.txt');

% Choose the image points
[H_affine, I_affine] = affine_rectification(I,line_info_file,614,159,645,541);

%% Compute the homography that metricly rectifies the image

[H_metric, I_metric] = metric_rectification(I_affine,H_affine,line_info_file,614,159,645,541);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)

I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% Indices of pairs of points per line
i = 493;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 48;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 186;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 508;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';
i = 424;
p9 = [A(i,1) A(i,2) 1]';
p10 = [A(i,3) A(i,4) 1]';
i = 712;
p11 = [A(i,1) A(i,2) 1]';
p12 = [A(i,3) A(i,4) 1]';
i = 240;
p13 = [A(i,1) A(i,2) 1]';
p14 = [A(i,3) A(i,4) 1]';
i = 565;
p15 = [A(i,1) A(i,2) 1]';
p16 = [A(i,3) A(i,4) 1]';
i = 357;
p17 = [A(i,1) A(i,2) 1]';
p18 = [A(i,3) A(i,4) 1]';
i = 119;
p19 = [A(i,1) A(i,2) 1]';
p20 = [A(i,3) A(i,4) 1]';

% Compute the lines that pass through the different pairs of points
l1 = cross(p1, p2);
m1 = cross(p3, p4);
l2 = cross(p5, p6);
m2 = cross(p7, p8);
l3 = cross(p9, p10);
m3 = cross(p11, p12);
l4 = cross(p13, p14);
m4 = cross(p15, p16);
l5 = cross(p17, p18);
m5 = cross(p19, p20);

[H_rect, I_rect] = single_step_rectification(I,l1,m1,l2,m2,l3,m3,l4,m4,l5,m5);
