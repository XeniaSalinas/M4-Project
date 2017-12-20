function [H_affine, I_affine] = affine_rectification(I, info_lines_file,...
pt_idx1,pt_idx2, pt_idx3, pt_idx4)

A = info_lines_file;
% Indices of pairs of points per line
i = pt_idx1;
p1 = [A(i,1) A(i,2) 1]';
p2 = [A(i,3) A(i,4) 1]';
i = pt_idx2;
p3 = [A(i,1) A(i,2) 1]';
p4 = [A(i,3) A(i,4) 1]';
i = pt_idx3;
p5 = [A(i,1) A(i,2) 1]';
p6 = [A(i,3) A(i,4) 1]';
i = pt_idx4;
p7 = [A(i,1) A(i,2) 1]';
p8 = [A(i,3) A(i,4) 1]';

% Compute the lines l1, l2, l3, l4, that pass through the different pairs of points
l1 = cross(p1, p2);
l2 = cross(p3, p4);
l3 = cross(p5, p6);
l4 = cross(p7, p8);

% Show the chosen lines in the image
figure, imshow(I), title('Original image');
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'g');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'b');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'r');

%Compute the homography that affinely rectifies the image

% Compute vanishing points
% cross-product between line 1 and 2
vp_1 = cross(l1, l2);
% cross-product between line 3 and 4
vp_2 = cross(l3, l4);

% Compute vanishing line
vline = cross(vp_1, vp_2);
vline = vline / vline(3);

% Construct affine recitification matrix
H_affine = eye(3);
H_affine(3,:) = vline;

% Apply affine rectification
I_affine = apply_H(I, H_affine);

% Compute the transformed lines lr1, lr2, lr3, lr4 (NOTE: l'=H-T *l)
% affine_rect_H_lines = transpose(inv(H_affine));
lr1 = inv(H_affine)' * l1;
lr2 = inv(H_affine)' * l2;
lr3 = inv(H_affine)' * l3;
lr4 = inv(H_affine)' * l4;

% show the transformed lines in the transformed image
figure, imshow(uint8(I_affine)), title('Image affinitely rectified');
hold on;
t=1:0.1:max(size(I_affine));
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
angleout = abs((a1>pi/2)*pi-a1)*180/pi;
disp(['Angle between two orthogonal lines before the affine recification: ', num2str(angleout)])

%Angle between two orthogonal lines after the affine recification
a1_affine = mod(atan2( det([lr1_e;lr3_e;]) , dot(lr1_e,lr3_e) ), 2*pi );
angleout_affine = abs((a1_affine>pi/2)*pi-a1_affine)*180/pi;
disp(['Angle between two orthogonal lines after the affine recification: ', num2str(angleout_affine)])

end