close all, clear
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
s = 0.8;                  % scale factor
angle = 30;               % rotation angle
t = [-1, 4];              % translation 
H = [s*cosd(angle) -s*sind(angle) t(1); ...
    s*sind(angle) s*cosd(angle) t(2); ...
    0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));


%% 1.2. Affinities

% Generate a matrix H which produces an affine transformation
A = [0.6 0.1;
     0.5 0.8];
T = [3; -4];
H = [A T; 0 0 1];

I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

% Decompose the affinity in four transformations: two
% rotations, a scale, and a translation
[U,D,V] = svd(A);
R_theta = U * transpose(V);
R_phi = transpose(V);
S = D(1:2,1:2);

% Verify that the product of the four previous transformations
% produces the same matrix H as above
tolerance = 1e-12;
A2 = R_theta * transpose(R_phi) * S * R_phi;
H2 = [A2 T; 0 0 1];
if abs(sum(sum(H - H2))) > tolerance
    error('H is no equal to its decomposition');
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
% I3 = apply_H(I, [zeros(2,2) T; 0 0 1]);
I3 = apply_H(I,[R_phi [0; 0]; 0 0 1]);
I3 = apply_H(I3,[S [0; 0]; 0 0 1]);
I3 = apply_H(I3,[transpose(R_phi) [0; 0]; 0 0 1]);
%Crop the last tranformation
I3 = apply_H(I3,[R_theta [0; 0]; 0 0 1]);
%Resize the transformation because the image sizes mismatch
% I3 = imresize(I3,[ size(I2,1)  size(I2,2)],'bilinear');
% % figure; imshow(uint8(I3));
% t2 = isequal(I2,I3)
% t2 = sum(sum((I2-I3)))

%Show the error
% figure; imshow(uint8(I2-I3)); 

if abs(sum(sum(I2 - I3))) > tolerance
    error('I3 is no equal to I2');
end

%% 1.3 Projective transformations (homographies)

% Generate a matrix H which produces a projective transformation
A = [0.9 0.12;
     0.05 1.1];
t = [3; -4];
v = [0.0005 0.001];
H = [A t; v 1];
I2 = apply_H(I, H);
figure; imshow(I); figure; imshow(uint8(I2));

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
I2 = apply_H(I, affine_rect_H);
figure; imshow(uint8(I2));

% Compute the transformed lines lr1, lr2, lr3, lr4 (NOTE: l'=H-T *l)
affine_rect_H_lines = transpose(inv(affine_rect_H));
lr1 = affine_rect_H_lines * l1;
lr2 = affine_rect_H_lines * l2;
lr3 = affine_rect_H_lines * l3;
lr4 = affine_rect_H_lines * l4;

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'g');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'b');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'r');

% To evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
l1_2 = [l1(1)/l1(3), l1(2)/l1(3)];
l2_2 = [l2(1)/l2(3), l2(2)/l2(3)];
l3_2 = [l3(1)/l3(3), l3(2)/l3(3)];
l4_2 = [l4(1)/l4(3), l4(2)/l4(3)];

lr1_2 = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
lr2_2 = [lr2(1)/lr2(3), lr2(2)/lr2(3)];
lr3_2 = [lr3(1)/lr3(3), lr3(2)/lr3(3)];
lr4_2 = [lr4(1)/lr4(3), lr4(2)/lr4(3)];

a1 = mod(atan2( det([l1_2;l2_2;]) , dot(l1_2,l2_2) ), 2*pi );
angleout = abs((a1>pi/2)*pi-a1);
a1_transf = mod(atan2( det([lr1_2;lr2_2;]) , dot(lr1_2,lr2_2) ), 2*pi );
angleout_transf = abs((a1_transf>pi/2)*pi-a1_transf);
angle_dif = (angleout - angleout_transf) * 180/pi

a2 = mod(atan2( det([l4_2;l3_2;]) , dot(l4_2,l3_2) ), 2*pi );
angleout = abs((a2>pi/2)*pi-a1);
a2_transf = mod(atan2( det([lr4_2;lr3_2;]) , dot(lr4_2,lr3_2) ), 2*pi );
angleout_transf = abs((a2_transf>pi/2)*pi-a2_transf);
angle_dif = (angleout - angleout_transf) * 180/pi


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

% Solve the system of equation provided by the orthogonal lines
A = [l(1)*m(1) l(1)*m(2)+l(1)*m(2)];
B = - l(2)*m(2);
s = linsolve(A,B);
S = [s(1) s(2); s(2) 1];

% Decompose the S matrix
[U,D,V] = svd(S);
A = U*sqrt(D(1:2,1:2))*transpose(U);
T = [0; 0];
H_2 = [A T; 0 0 1];

%Restore the image
H_metric = inv(H_2);
I_metric = apply_H(I_affine, H_metric);
figure; imshow(uint8(I_metric));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



