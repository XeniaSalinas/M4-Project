function v_t = compute_row_absolute_conic(H, i, j)    
%COMPUTE_ROW_ABSOLUTE_CONIC Computes Vt(i,j) given homography H and a
%position i, j. Used to obtain the image of the absolute conic w
h1 = H(1,:);
h2 = H(2,:);
h3 = H(3,:);
v_t = [h1(i)*h1(j), h1(i)*h2(j)+h2(i)*h1(j), h1(i)*h3(j)+h3(i)*h1(j), ...
       h2(i)*h2(j), h2(i)*h3(j)+h3(i)*h2(j), h3(i)*h3(j)];
end  

