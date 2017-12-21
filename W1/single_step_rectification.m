function [H_rect, I_rect] = single_step_rectification(I,l1,m1,l2,m2,l3,m3,l4,m4,l5,m5)
%SINGLE_STEP_RECTIFICATION Metric Rectification in a single step

% Show the chosen lines in the image
figure, imshow(I), title('Original image');
hold on;
t=1:0.1:1000;
plot(t, -(l1(1)*t + l1(3)) / l1(2), 'y');
plot(t, -(m1(1)*t + m1(3)) / m1(2), 'y');
plot(t, -(l2(1)*t + l2(3)) / l2(2), 'b');
plot(t, -(m2(1)*t + m2(3)) / m2(2), 'b');
plot(t, -(l3(1)*t + l3(3)) / l3(2), 'r');
plot(t, -(m3(1)*t + m3(3)) / m3(2), 'r');
plot(t, -(l4(1)*t + l4(3)) / l4(2), 'g');
plot(t, -(m4(1)*t + m4(3)) / m4(2), 'g');
plot(t, -(l5(1)*t + l5(3)) / l5(2), 'c');
plot(t, -(m5(1)*t + m5(3)) / m5(2), 'c');

% Solve the system of equation provided by the orthogonal lines
B = [l1(1)*m1(1), (l1(1)*m1(2)+l1(2)*m1(1))/2, l1(2)*m1(2), (l1(1)*m1(3)+l1(3)*m1(1))/2, (l1(2)*m1(3)+l1(3)*m1(2))/2, l1(3)*m1(3);
     l2(1)*m2(1), (l2(1)*m2(2)+l2(2)*m2(1))/2, l2(2)*m2(2), (l2(1)*m2(3)+l2(3)*m2(1))/2, (l2(2)*m2(3)+l2(3)*m2(2))/2, l2(3)*m2(3);
     l3(1)*m3(1), (l3(1)*m3(2)+l3(2)*m3(1))/2, l3(2)*m3(2), (l3(1)*m3(3)+l3(3)*m3(1))/2, (l3(2)*m3(3)+l3(3)*m3(2))/2, l3(3)*m3(3);
     l4(1)*m4(1), (l4(1)*m4(2)+l4(2)*m4(1))/2, l4(2)*m4(2), (l4(1)*m4(3)+l4(3)*m4(1))/2, (l4(2)*m4(3)+l4(3)*m4(2))/2, l4(3)*m4(3);
     l5(1)*m5(1), (l5(1)*m5(2)+l5(2)*m5(1))/2, l5(2)*m5(2), (l5(1)*m5(3)+l5(3)*m5(1))/2, (l5(2)*m5(3)+l5(3)*m5(2))/2, l5(3)*m5(3)];
s = null(B); % Null vector of B.
C = [s(1), s(2)/2, s(4)/2;
     s(2)/2, s(3), s(5)/2;
     s(4)/2, s(5)/2, s(6)];

[U,D,~] = svd(C);
H_rect = inv(U*D);

%Restore the image
I_rect = apply_H(I, H_rect);

if (I_rect(1,1,1) == 0 && I_rect(1,1,2) == 0 && I_rect(1,1,3) == 0)
    offset = find(I_rect(1,:,1)~=0,1);
    H_rect(1,3) = offset;
end

% Compute the transformed lines lr, mr (NOTE: l'=H-T *l)
lr1 = inv(H_rect)' * l1;
mr1 = inv(H_rect)' * m1;
lr2 = inv(H_rect)' * l2;
mr2 = inv(H_rect)' * m2;

% show the transformed lines in the transformed image
figure, imshow(uint8(I_rect));
title('Image metricly rectified');
hold on;
t=1:0.1:max(size(I_rect));
plot(t, -(lr1(1)*t + lr1(3)) / lr1(2), 'y');
plot(t, -(mr1(1)*t + mr1(3)) / mr1(2), 'g');
plot(t, -(lr2(1)*t + lr2(3)) / lr2(2), 'b');
plot(t, -(mr2(1)*t + mr2(3)) / mr2(2), 'r');

end

