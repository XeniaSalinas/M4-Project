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


%% 1.1. Similarities
I=imread('Data/0005_s.png'); % we have to be in the proper folder

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
t = [3; -4];
H = [A t; 0 0 1];

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
H2 = [A2 t; 0 0 1];
if abs(sum(sum(H - H2))) > tolerance
    error('H is no equal to its decomposition');
end

% ToDo: verify that the proper sequence of the four previous
% transformations over the image I produces the same image I2 as before
I3 = apply_H(I, [zeros(2,2) T; 0 0 1]);
I3 = apply_H(I,[R_phi [0; 0]; 0 0 1]);
I3 = apply_H(I3,[S [0; 0]; 0 0 1]);
I3 = apply_H(I3,[transpose(R_phi) [0; 0]; 0 0 1]);
I3 = apply_H(I3,[R_theta [0; 0]; 0 0 1]);
t2 = isequal(I2,I3)
t2 = sum(sum((I2-I3)))

figure; imshow(uint8(I2)); 

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


% choose the image points
I=imread('Data/0000_s.png');
A = load('Data/0000_s_info_lines.txt');

% indices of lines
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
l3 = cross(p5,p6);
l4 = cross(p7, p8);

% show the chosen lines in the image
figure;imshow(I);
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'r');

% ToDo: compute the homography that affinely rectifies the image
%%
I2 = apply_H(I, H);
%figure; imshow(uint8(I2));

% ToDo: compute the transformed lines lr1, lr2, lr3, lr4 (NOTE: l'=H-T *l)
lr1 = inv(H).' * l1;
lr2 = inv(H).' * l2;
lr3 = inv(H).' * l3;
lr4 = inv(H).' * l4;

% show the transformed lines in the transformed image
figure;imshow(uint8(I2));
hold on;
t=1:0.1:1000;
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'g');
plot(t, -(lr3(1)*t + lr3(3)) / lr3(2), 'b');
plot(t, -(lr4(1)*t + lr4(3)) / lr4(2), 'r');

% ToDo: to evaluate the results, compute the angle between the different pair 
% of lines before and after the image transformation
% l1_2 = [l1(1)/l1(3), l1(2)/l1(3)];
% l2_2 = [l2(1)/l2(3), l2(2)/l2(3)];
% l3_2 = [l3(1)/l3(3), l3(2)/l3(3)];
% l4_2 = [l4(1)/l4(3), l4(2)/l4(3)];
% 
% lr1_2 = [lr1(1)/lr1(3), lr1(2)/lr1(3)];
% lr2_2 = [lr2(1)/lr2(3), lr2(2)/lr2(3)];
% lr3_2 = [lr3(1)/lr3(3), lr3(2)/lr3(3)];
% lr4_2 = [lr4(1)/lr4(3), lr4(2)/lr4(3)];
% 
% a1 = mod(atan2( det([l1_2;l2_2;]) , dot(l1_2,l2_2) ), 2*pi );
% angleout = abs((a1>pi/2)*pi-a1);
% a1_transf = mod(atan2( det([lr1_2;lr2_2;]) , dot(lr1_2,lr2_2) ), 2*pi );
% angleout_transf = abs((a1_transf>pi/2)*pi-a1_transf);
% angle_dif = (angleout - angleout_transf) * 180/pi
% 
% a2 = mod(atan2( det([l4_2;l3_2;]) , dot(l4_2,l3_2) ), 2*pi );
% angleout = abs((a2>pi/2)*pi-a1);
% a2_transf = mod(atan2( det([lr4_2;lr3_2;]) , dot(lr4_2,lr3_2) ), 2*pi );
% angleout_transf = abs((a2_transf>pi/2)*pi-a2_transf);
% angle_dif = (angleout - angleout_transf) * 180/pi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Metric Rectification

%% 3.1 Metric rectification after the affine rectification (stratified solution)

% ToDo: Metric rectification (after the affine rectification) using two non-parallel orthogonal line pairs
%       As evaluation method you can display the images (before and after
%       the metric rectification) with the chosen lines printed on it.
%       Compute also the angles between the pair of lines before and after
%       rectification.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Affine and Metric Rectification of the left facade of image 0001

% ToDo: Write the code that rectifies the left facade of image 0001 with
%       the stratified method (affine + metric). 
%       Crop the initial image so that only the left facade is visible.
%       Show the (properly) transformed lines that use in every step.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. OPTIONAL: Metric Rectification in a single step
% Use 5 pairs of orthogonal lines (pages 55-57, Hartley-Zisserman book)



