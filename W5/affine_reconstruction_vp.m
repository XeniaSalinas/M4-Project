function [Hp] = affine_reconstruction_vp(vanishing_points_1, vanishing_points_2, P1, P2, w, h)
%AFFINE_RECONSTRUCTION Return affine reconstruction matrix Hp given
% vanishing points from 2 images, camera matrices and width and height of
% the images.

v1 = vanishing_points_1{1};
v2 = vanishing_points_1{2};
v3 = vanishing_points_1{3};

v1p = vanishing_points_2{1};
v2p = vanishing_points_2{2};
v3p = vanishing_points_2{3};

A = [...
    triangulate(euclid(v1), euclid(v1p), P1, P2, [w,h])';
    triangulate(euclid(v2), euclid(v2p), P1, P2, [w,h])';
    triangulate(euclid(v3), euclid(v3p), P1, P2, [w,h])';...
    ];

p = null(A);
p = p/p(end);
Hp = [eye(3) zeros(3,1);
      transpose(p)];

end

