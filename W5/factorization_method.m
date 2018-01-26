function [Pproj, Xproj] = factorization_method(x)

    tol_d = 0.1;
    d = inf;
    
    num_cameras = size(x,2);
    num_points = size(x{1},2);
    
    % Initialize lambda
    lambda = ones(num_cameras, num_points);
    
    % 1. Normalize the image coordinates, by applying transformations Ti
    x_hat = cell(1,num_cameras);
    T = cell(1,num_cameras);
    for i = 1:num_cameras
        [x_hat{i}, T{i}] = normalise2dpts(x{i});
    end
    
    % 2. Estimate the fundamental matrices and epipoles with the method of
    % [Hartley, 1995]
    F = cell(1, num_cameras);
    e = cell(1, num_cameras);
    for i = 1:num_cameras
        F{i} = fundamental_matrix(x{i}, x{1});
        [U, D, V] = svd(F{i});
        e{i} = V(:,3) / V(3,3); 
    end
    
    % 3. Determine the scale factors lambda_ip using equation (3) of [Sturm and Triggs, 1996]
    for j=1:size(x{1},2)
        num = x{1}(:, j)' * F{1} * cross(e{1}, x{1}(:,j));
        den = norm(cross(e{1}, x{1}(:,j))) .^ 2 * lambda(1, j);
        lambda(1,j) = num / den;
    end
    for j=1:size(x{2},2)
        num = x{1}(:, j)' * F{2} * cross(e{2}, x{2}(:,j));
        den = norm(cross(e{2}, x{2}(:,j))) .^ 2 * lambda(1, j);
        lambda(2,j) = num / den;
    end
    
    % 4. Build the rescaled measurement matrix W
    
    % 5. Balance W by column-wise and "triplet-o-rows"-wise scalar
    % multiplications
    
    % 6. compute the SVD of the balanced matrix W
    
    % 7. From the SVD, recover projective motion and shape
    
    % 8. Adapt projective motion, to account for the normalization
    % transformations Ti os step 1
    
    true = 1;
    while(true)

    end

end
