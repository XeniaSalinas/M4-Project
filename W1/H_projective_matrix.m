function H = H_projective_matrix(affine_transformation, translation, translation_vector)
% Compute an affine matrix
H = [affine_transformation translation;translation_vector 1];
end