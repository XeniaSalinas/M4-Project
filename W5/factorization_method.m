function [Pproj, Xproj, reproj_err] = factorization_method(x, initialization)
% TODO: Extend this function to use more than 2 cameras. Now it assumes 2
% cameras are passed.

tol_d = 0.1;
d = inf;

num_cameras = size(x,2);
if num_cameras ~= 2
    error('Projective reconstruction (factorization) not implemented for more than 2 cameras');
end
num_points = size(x{1},2);

% Initialize lambda
temp_lambda = ones(num_cameras, num_points);

% 1. Normalize the image coordinates, by applying transformations Ti
x_hat = cell(1,num_cameras);
T = cell(1,num_cameras);
for i = 1:num_cameras
    [x_hat{i}, T{i}] = normalise2dpts(x{i});
end
% Sturm initialization of lambdas
if isequal(initialization, 'sturm')
    % 2. Estimate the fundamental matrix and epipoles with the method of
    % [Hartley, 1995]
    F = cell(1, num_cameras);
    e = cell(1, num_cameras);
    for i = 2:num_cameras
        F{i} = ransac_fundamental_matrix(x_hat{i}, x_hat{i-1}, 2.0);
        [~, ~, V] = svd(F{i});
        e{i} = V(:,3) / V(3,3); 
    end

    % 3. Determine the scale factors lambda_ip using equation (3) of [Sturm and Triggs, 1996]
    for i=2:num_cameras
       e_rep = repmat(e{i}, 1, num_points);
       xprod = cross(e_rep, x_hat{i});
       num = diag((x_hat{i-1}' * F{i}) * xprod)';
       den = sum(xprod.^2, 1);
       temp_lambda(i,:) = (num ./ den) .* temp_lambda(i-1, :);
    end
 end
count = 1;

while(1)
    temp_diff = Inf;
    norm_row = false;
    % 4. Find iteratively correct labmdas to build the rescaled measurement matrix W
    while(1)
        norm_row =~norm_row;
        diff = temp_diff;
        lambda = temp_lambda;
        
        % 5. Balance W by column-wise and "triplet-o-rows"-wise scalar
        % multiplications
        if norm_row
            temp_lambda(1,:) = temp_lambda(1,:) ./ norm(temp_lambda(1,:));
            temp_lambda(2,:) = temp_lambda(2,:) ./ norm(temp_lambda(2,:));
        else
            for col = 1:size(x{1},2)
                temp_lambda(:,col) = temp_lambda(:,col) / norm(temp_lambda(:,col));
            end
        end
        
        temp_diff = (lambda - temp_lambda).^2;
        temp_diff = sum(temp_diff(:));
        temp_diff = sqrt(temp_diff);

        if temp_diff == 0 || ((abs(temp_diff - diff)/temp_diff) < tol_d)
            lambda = temp_lambda;
            break;
        end

    end
    % Build the rescaled measurement matrix W
    W = zeros(3*2, size(x{1},2));
    W(1,:) = lambda(1,:) .* x_hat{1}(1,:);
    W(2,:) = lambda(1,:) .* x_hat{1}(2,:);
    W(3,:) = lambda(1,:) .* x_hat{1}(3,:);
    W(4,:) = lambda(2,:) .* x_hat{2}(1,:);
    W(5,:) = lambda(2,:) .* x_hat{2}(2,:);
    W(6,:) = lambda(2,:) .* x_hat{2}(3,:);
    % 6. compute the SVD of the balanced matrix W
    [U,D,V]=svd(W);
    Xproj=V(:,1:4)';
    P_motion=U*D(:,1:4);

    d_old = d;
    d = 0;

    %Compute distance from the real point to the 3D reconstructed point
    Px1 = P_motion(1:3,:) * Xproj;
    Px1 = Px1./Px1(end,:);
    x_1 = x_hat{1};
    Px2 = P_motion(4:6,:) * Xproj;
    Px2 = Px2./Px2(end,:);
    x_2 = x_hat{2};
    for j=1:size(x{1},2)
         d = d + sum((x_1(:,j) - Px1(:,j)).^2)+ sum((x_2(:,j) - Px2(:,j)).^2);
    end

    if ((abs(d - d_old)/d) < tol_d)
        count
        break;
    else
        %Update lambdas
        temp = P_motion*Xproj;
        temp_lambda(1,:) = temp(3,:);
        temp_lambda(2,:) = temp(6,:);
    end
    count=count+1;
end

% 8. Adapt projective motion, to account for the normalization
% transformations Ti of step 1
Pproj(1:3,:) = T{1}\ P_motion(1:3,:);
Pproj(4:6,:) = T{2}\ P_motion(4:6,:);
reproj_err = d; % Reprojection error

end
