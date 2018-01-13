function F = fundamental_matrix(x1_test, x2_test)

% Normalize the input points and get H and H'
[x1_norm, H1] = normalise2dpts(x1_test);
[x2_norm, H2] = normalise2dpts(x2_test);

x1 = x1_norm(1,:)';
y1 = x1_norm(2,:)';

x2 = x2_norm(1,:)';
y2 = x2_norm(2,:)';
%Total points (8 correspondences)
N = size(x1_norm, 2);

%W = zeros(N,9);
W = [x1.*x2 y1.*x2 x2 x1.*y2 y1.*y2 y2 x1 y1 ones(N, 1)];

% Singular Value Decomposition of W
[~, ~, V] = svd(W);

% F is the last column of V
f = V(:,end);
% Reshape as a matrix
F = [f(1) f(2) f(3); f(4) f(5) f(6); f(7) f(8) f(9)];

% Singular Value Decomposition of F
[U, D, V] = svd(F);

% Delete the minimum singular value of D to make it rank 2
D_hat = D;
D_hat(3,3) = 0;

% Recompute F of rank 2
F = U * D_hat * V';

% Denormalize F using the 3x3 tranformation matrices of x1 and x2
F = H2' * F * H1; 

end