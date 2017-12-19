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
s = 1;                  % scale factor
angle = 30;               % rotation angle
t = [-1, 4];              % translation 
H = [s*cosd(angle) -s*sind(angle) t(1); ...
    s*sind(angle) s*cosd(angle) t(2); ...
    0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% Generate a matrix H which produces an affine transformation
A = [1 -0.1;
     0.5 0.8];
T = [20; 30];
H = [A T; 0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

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
H_decomposed = [A2 T; 0 0 1];
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
figure; imshow(uint8(I3));

% Compare the image transformed using H and the image transformed using the
% decomposition of H.
[rows, cols, c] = size(I2);
I3_resize = imresize(I3, [rows, cols]);
error_im = uint8(I2 - I3_resize);
mae = mean2(error_im);  % mean absolute difference
fprintf('Mean average error: %.3f\n', mae);
figure; imshow(error_im); 

%% 1.3 Projective transformations (homographies)

% Generate a matrix H which produces a projective transformation
A = [0.9 0.12;
     0.05 1.1];
t = [3; -4];
v = [0.0005 0.001];
H_projective = [A t; v 1];
I_projective = apply_H(I, H_projective);
figure; imshow(I); figure; imshow(uint8(I_projective));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Affine Rectification


% Choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% Indices of pairs of points per line
i = 424;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 240;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 712;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 565;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% Compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1, p2);
l2 = cross(p3, p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);

% Show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'r');

%% Compute the homography that affinely rectifies the image

% Compute vanishing points
% cross-product between line 424 and 240
vp_1 = cross(l1, l2);
% cross-product between line 712 and 565
vp_2 = cross(l3, l4);

% Compute vanishing line
vline = cross(vp_1, vp_2);
vline = vline / vline(3);

% Construct affine recitification matrix
affine_rect_H = eye(3);
affine_rect_H(3,:) = vline;

% Apply affine rectification
I_affine = apply_H(I, affine_rect_H);
figure; imshow(uint8(I_affine));

% Compute the transformed lines lr1, lr2, lr3, lr4 (NOTE: l'=H-T *l)
affine_rect_H_lines = transpose(inv(affine_rect_H));
lr1 = affine_rect_H_lines * l1;
lr2 = affine_rect_H_lines * l2;
lr3 = affine_rect_H_lines * l3;
lr4 = affine_rect_H_lines * l4;

% show the transformed lines in the transformed image
figure;imshow(uint8(I_affine));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'g');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'b');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'r');

% To evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation

%Euclidean representation of the lines before the transformation
l1_e = [l1(1)/l1(3), l1(2)/l1(3)];
l3_e = [l3(1)/l3(3), l3(2)/l3(3)];

%Euclidean representation of the lines after the transformation
lr1_e = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
lr3_e = [lr3(1)/lr3(3), lr3(2)/lr3(3)];


%Angle between two orthogonal lines before the affine recification
a1 = mod(atan2( det([l1_e;l3_e;]) , dot(l1_e,l3_e) ), 2*pi );
angleout = abs((a1>pi/2)*pi-a1)*180/pi

%Angle between two orthogonal lines after the affine recification
a1_affine = mod(atan2( det([lr1_e;lr3_e;]) , dot(lr1_e,lr3_e) ), 2*pi );
angleout_affine = abs((a1_affine>pi/2)*pi-a1_affine)*180/pi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%Choose two orthogonal lines
l = lr1;
m = lr3;

% Get the points in the corners of the window:
x1 = cross(lr1, lr3);
x2 = cross(lr2, lr4);
x3 = cross(lr1, lr4);
x4 = cross(lr2, lr3);
% With these points, compute the diagonal lines:
v1 = cross(x1, x2);
v2 = cross(x3, x4);

% Solve the system of equation provided by the orthogonal lines
B = [l(1)*m(1), l(1)*m(2)+l(2)*m(1), l(2)*m(2);
    v1(1)*v2(1), v1(1)*v2(2)+v1(2)*v2(1), v1(2)*v2(2)];
s = null(B); % Null vector of B.
S = [s(1), s(2); s(2), s(3)];

% Compute the upper triangular matrix using the Cholesky factorization:
K = inv(chol(inv(S)));

T = [0; 0];
H_2 = [K T; 0 0 1];
H_metric = inv(H_2);

%Restore the image
I_metric = apply_H(I_affine, H_metric);

% Compute the transformed lines lr, mr (NOTE: l'=H-T *l)
metric_rect_H_lines = transpose(inv(H_metric));
lr = metric_rect_H_lines * l;
mr = metric_rect_H_lines * m;

% show the transformed lines in the transformed image
figure; imshow(uint8(I_metric));
hold on;
t=1:0.1:1000;
plot(t, -(lr(1)*t + lr(3)) / lr(2), 'y');
plot(t, -(mr(1)*t + mr(3)) / mr(2), 'g');

%Euclidean coordinates of the lines before the metric rectification
l_e = [l(1)/l(3), l(2)/l(3)];
m_e = [m(1)/m(3), m(2)/m(3)];

%Euclidean coordinates of the lines after the metric rectification
lr_e = [lr(1)/lr(3), lr(2)/lr(3)];
mr_e = [mr(1)/mr(3), mr(2)/mr(3)];

%Angle between two orthogonal lines before the metric recification
a2 = mod(atan2( det([m_e;l_e;]) , dot(m_e,l_e) ), 2*pi );
angleout = abs((a2>pi/2)*pi-a1)* 180/pi

%Angle between two orthogonal lines after the metric recification
a2_metric = mod(atan2( det([mr_e;lr_e;]) , dot(mr_e,lr_e) ), 2*pi );
angleout_metric = abs((a2_metric>pi/2)*pi-a2_metric)* 180/pi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

I=imread('Data/0001_s.png');
I = imcrop(I,[1 1 533 614]);

A = load('Data/0001_s_info_lines.txt');

% Indices of pairs of points per line
i = 614;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = 159;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = 645;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = 541;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% Compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1, p2);
l2 = cross(p3, p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);

% Show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'r');

%% Compute the homography that affinely rectifies the image

% Compute vanishing points
% cross-product between line 424 and 240
vp_1 = cross(l1, l2);
% cross-product between line 712 and 565
vp_2 = cross(l3, l4);

% Compute vanishing line
vline = cross(vp_1, vp_2);
vline = vline / vline(3);

% Construct affine recitification matrix
affine_rect_H = eye(3);
affine_rect_H(3,:) = vline;

% Apply affine rectification
I_affine = apply_H(I, affine_rect_H);
figure; imshow(uint8(I_affine));

% Compute the transformed lines lr1, lr2, lr3, lr4 (NOTE: l'=H-T *l)
affine_rect_H_lines = transpose(inv(affine_rect_H));
lr1 = affine_rect_H_lines * l1;
lr2 = affine_rect_H_lines * l2;
lr3 = affine_rect_H_lines * l3;
lr4 = affine_rect_H_lines * l4;

% show the transformed lines in the transformed image
figure;imshow(uint8(I_affine));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'g');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'b');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'r');

%% Compute the homography that metricly rectifies the image
%Choose two orthogonal lines
l = lr1;
m = lr3;

% Get the points in the corners of the window:
x1 = cross(lr1, lr3);
x2 = cross(lr2, lr4);
x3 = cross(lr1, lr4);
x4 = cross(lr2, lr3);
% With these points, compute the diagonal lines:
v1 = cross(x1, x2);
v2 = cross(x3, x4);

% Solve the system of equation provided by the orthogonal lines
B = [l(1)*m(1), l(1)*m(2)+l(2)*m(1), l(2)*m(2);
    v1(1)*v2(1), v1(1)*v2(2)+v1(2)*v2(1), v1(2)*v2(2)];
s = null(B); % Null vector of B.
S = [s(1), s(2); s(2), s(3)];

% Compute the upper triangular matrix using the Cholesky factorization:
K = inv(chol(inv(S)));

T = [0; 0];
H_2 = [K T; 0 0 1];
H_metric = inv(H_2);

%Restore the image
I_metric = apply_H(I_affine, H_metric);

% Compute the transformed lines lr, mr (NOTE: l'=H-T *l)
metric_rect_H_lines = transpose(inv(H_metric));
lr = metric_rect_H_lines * l;
mr = metric_rect_H_lines * m;

% show the transformed lines in the transformed image
figure; imshow(uint8(I_metric));
hold on;
t=1:0.1:3000;
plot(t, -(lr(1)*t + lr(3)) / lr(2) -1000, 'y');
plot(t, -(mr(1)*t + mr(3)) / mr(2) +1850, 'g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



