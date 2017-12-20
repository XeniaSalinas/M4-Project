function H = H_affine_matrix(affine_transformation, translation)
% Compute an affine matrix
H = [affine_transformation translation;0 0 1];
end