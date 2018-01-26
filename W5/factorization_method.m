function [Pproj, Xproj] = factorization_method(x)

    tol_d = 0.1;
    d = inf;
    
    num_cameras = size(x,2);
    num_points = size(x{1},2);
    
    % Initialize lambda
    lambda = ones(num_cameras, num_points);

    % Sturm initialization   
    F = cell(1, num_cameras);
    e = cell(1, num_cameras);
    for i = 1:num_cameras
        F{i} = fundamental_matrix(x{i}, x{1});
        [U, D, V] = svd(F{i});
        e{i} = V(:,3) / V(3,3); 
    end
    
    for c = 1:num_cameras
        for i=1:num_points
            lambda(1,i) = (x{1}(:, i)' * F{c} * cross(e{c}, x{c}(:,i))) / (norm(cross(e{c}, x{c}(:,i))).^ 2 * lambda(1, i));
        end
    end
    
    % Normalize input points
    x_hat = cell(1,num_cameras);
    T = cell(1,num_cameras);
    for i = 1:num_cameras
        [x_hat{i}, T{i}] = normalise2dpts(x{i});
    end
    
    true = 1;
    while(true)


    end

end
