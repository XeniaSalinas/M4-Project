function [ X ] = triangulate( x1, x2, P1, P2, imsize )
%TRIANGULATE Summary of this function goes here
%   Detailed explanation goes here

% %Homogeneous coordinates
% x1 = [x1; 1];
% x2 = [x2; 1];

%Normalization
H = [2/imsize(1) 0 -1;
    0 2/imsize(2) -1;
    0 0 1];

x1 = H*[x1' 1]';
x2 = H*[x2' 1]';
P1 = H*P1;
P2 = H*P2;

%Triangulation
A = [x1(1) * P1(3,:) - P1(1,:);
    x1(2) * P1(3,:) - P1(2,:);
    x2(1) * P2(3,:) - P2(1,:);
    x2(2) * P2(3,:) - P2(2,:)];

[U,D,V] = svd(A);

X = V(:,end);
X = X/X(end);
end

