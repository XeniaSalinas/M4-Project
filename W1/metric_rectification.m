function [H_metric, I_metric] = metric_rectification(I_affine, H_affine, info_lines_file,...
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

% Transform the lines:
l1_a = inv(H_affine)' * cross(p1, p2);
l2_a = inv(H_affine)' * cross(p3, p4);
l3_a = inv(H_affine)' * cross(p5, p6);
l4_a = inv(H_affine)' * cross(p7, p8);

%Choose two orthogonal lines
l = l1_a;
m = l3_a;

% Get the points in the corners of the window:
x1 = cross(l1_a, l3_a);
x2 = cross(l2_a, l4_a);
x3 = cross(l1_a, l4_a);
x4 = cross(l2_a, l3_a);

% With these points, compute the diagonal lines:
l2 = cross(x1, x2);
m2 = cross(x3, x4);
    
% Solve the system of equation provided by the orthogonal lines
B = [l(1)*m(1), l(1)*m(2)+l(2)*m(1), l(2)*m(2);
    l2(1)*m2(1), l2(1)*m2(2)+l2(2)*m2(1), l2(2)*m2(2)];
s = null(B); % Null vector of B.
S = [s(1), s(2); s(2), s(3)];

% Compute the upper triangular matrix using the Cholesky factorization:
K = inv(chol(inv(S))); %affine transformation
T = [0; 0]; %translation

H_2 = H_affine_matrix(K, T);
H_metric = inv(H_2);

%Restore the image
I_metric = apply_H(I_affine, H_metric);

if (I_metric(1,1,1) == 0 && I_metric(1,1,2) == 0 && I_metric(1,1,3) == 0)
    offset = find(I_metric(1,:,1)~=0,1);
    H_metric(1,3) = offset;
end

% Compute the transformed lines lr, mr (NOTE: l'=H-T *l)
lr = inv(H_metric)' * l;
mr = inv(H_metric)' * m;
lr2 = inv(H_metric)' * l2;
mr2 = inv(H_metric)' * m2;

% show the transformed lines in the transformed image
figure, imshow(uint8(I_metric));
title('Image metricly rectified');
hold on;
t=1:0.1:max(size(I_metric));
plot(t, -(lr(1)*t + lr(3)) / lr(2), 'y');
plot(t, -(mr(1)*t + mr(3)) / mr(2), 'g');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'b');
plot(t, -(mr2(1)*t + mr2(3)) / mr2(2), 'r');

%Euclidean coordinates of the lines before the metric rectification
l_e = [l(1)/l(3), l(2)/l(3)];
m_e = [m(1)/m(3), m(2)/m(3)];

%Euclidean coordinates of the lines after the metric rectification
lr_e = [lr(1)/lr(3), lr(2)/lr(3)];
mr_e = [mr(1)/mr(3), mr(2)/mr(3)];

%Angle between two orthogonal lines before the metric recification
a2 = mod(atan2( det([m_e;l_e;]) , dot(m_e,l_e) ), 2*pi );
angleout = abs((a2>pi/2)*pi-a2)* 180/pi;
disp(['Angle between two orthogonal lines before the metric recification: ', num2str(angleout)])

%Angle between two orthogonal lines after the metric recification
a2_metric = mod(atan2( det([mr_e;lr_e;]) , dot(mr_e,lr_e) ), 2*pi );
angleout_metric = abs((a2_metric>pi/2)*pi-a2_metric)* 180/pi;
disp(['Angle between two orthogonal lines after the metric recification: ', num2str(angleout_metric)])

end