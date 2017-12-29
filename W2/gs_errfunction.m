function [ Y_initial ] = gs_errfunction( P0, Xobs )
%GS_ERRFUNCTION

% get the homography
H = reshape(P0(1:9), [3,3]);

N = size(Xobs,1) / 2;

x = Xobs(1:N);
x = reshape(x, [2,size(x,1)/2]);

xp = Xobs(N+1:end);
xp = reshape(xp, [2,size(xp,1)/2]);

%Discard the Homography and get x_hat
x_hat = P0(10:end);
len_x_hat = size(x_hat,1)/2;
x_hat = reshape(x_hat, [2,len_x_hat]);
x_hat = [x_hat; ones(1,len_x_hat)];

% Calculate x_hat_p
x_hat_p = H * x_hat;

% Error between estimates and observations
error_x_x_hat = abs(x - euclid(x_hat));
error_xp_x_hat_p = abs(xp - euclid(x_hat_p));

Y_initial = [error_x_x_hat(:) error_xp_x_hat_p(:)];

end

