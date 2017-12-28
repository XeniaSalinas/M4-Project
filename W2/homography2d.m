function [ H ] = homography2d( x1, x2 )
%HOMOGRAPHY2D Summary of this function goes here
%   Detailed explanation goes here
% x1:  input points
% x2: output corresponding points 
% H: homography estimated with the correspondences x1 and x2
%First step: normalise data (as in DLT algorithm)
[x1, T] = DLT_normalization(x1);
[x2, T_out] = DLT_normalization(x2);

%Second step: 
n = size(x1, 2);
x = x2(1, :); 
y = x2(2,:); 
X = x1(1,:); 
Y = x1(2,:);
row_zeros = zeros(3, n);
row_values = -[X; Y; ones(1,n)];
hx = [row_values; row_zeros; x.*X; x.*Y; x];
hy = [row_zeros; row_values; y.*X; y.*Y; y];
h = [hx hy];
[~, ~, V] = svd(h');
v = V(:,9);
H = transpose(reshape(v(:,1), [3, 3]));
H = inv(T_out)*H*T;

end