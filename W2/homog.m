function [ x_homog ] = homog( x_euclidean )
x_homog = [x_euclidean; ones(1, size(x_euclidean,2))];
end