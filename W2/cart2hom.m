function x_h = cart2hom(x_c)
%CART2HOM Cartesian (inhomogeneous) coordinates to homogeneous coordinates

% Expects N points expressed in cartesian coordinates as a NxD matrix,
% where D is the number of dimensions of the cartesian space (tipically 2
% or 3).
[N, D] = size(x_c);
x_h = [x_c, ones(N, 1)];

end

