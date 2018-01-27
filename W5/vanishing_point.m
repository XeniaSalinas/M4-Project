function [ v ] = vanishing_point( xo1, xf1, xo2, xf2 )
%VANISHING_POINT Summary of this function goes here
%   Detailed explanation goes here
l1 = cross(xo1,xf1);
l2 = cross(xo2, xf2);
l1 = l1/l1(end);
l2 = l2/l2(end);
v = cross(l1,l2);

end
