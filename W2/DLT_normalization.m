function [ xn, T ] = DLT_normalization( x )
%DLT_NORMALIZATION Summary of this function goes here
%   Detailed explanation goes here
    
    centroid=mean(x(1:end-1,:)')';
    x_cent(1,:) =  x(1,:)-centroid(1);
    x_cent(2,:) = x(2,:)-centroid(2);
    
    scale = sqrt(2) / mean( sqrt(sum(x_cent(1:end-1,:).^2)) );
    T = [scale   0   -scale*centroid(1)
         0     scale -scale*centroid(2)
         0       0      1      ];
    xn = T*x;

end

