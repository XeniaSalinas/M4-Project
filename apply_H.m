function I2 = apply_H(I, H)
%APPLY_H Transform the input image with the given homography

%Compute the size of the transformed image
[m,n,c] = size(I);
lu_p = [1;1;1];
ru_p = [m;1;1];
ld_p = [1;n;1];
rd_p = [m;n;1];

lu2_p = H * lu_p;
ru2_p= H * ru_p;
ld2_p = H * ld_p;
rd2_p = H * rd_p;

lu2_e = [lu2_p(1)/lu2_p(3), lu2_p(2)/lu2_p(3)];
ru2_e = [ru2_p(1)/ru2_p(3), ru2_p(2)/ru2_p(3)];
ld2_e = [ld2_p(1)/ld2_p(3), ld2_p(2)/ld2_p(3)];
rd2_e = [rd2_p(1)/rd2_p(3), rd2_p(2)/rd2_p(3)];

m2 = ceil(abs(max([lu2_e(1), ru2_e(1), ld2_e(1), rd2_e(1)]) - min([lu2_e(1), ru2_e(1), ld2_e(1), rd2_e(1)]))) + 1;
n2 = ceil(abs(max([lu2_e(2), ru2_e(2), ld2_e(2), rd2_e(2)]) - min([lu2_e(2), ru2_e(2), ld2_e(2), rd2_e(2)]))) + 1;

%Compute the position rectification vector
origin_m = min([lu2_e(1), ru2_e(1), ld2_e(1), rd2_e(1)]);
origin_n = min([lu2_e(2), ru2_e(2), ld2_e(2), rd2_e(2)]);

%Compute the homography transfromation
I2 = zeros(m2,n2,c);

for i = 1:m2
    for j = 1:n2
        x2_e = [i; j] + [origin_m; origin_n];
        x2_p = H\[x2_e;1]; %A\b for inv(a)*b
        x_e = round([x2_p(1)/x2_p(3), x2_p(2)/x2_p(3)]); %%TODO: Use interp2
        if (x_e(1) > 1 && x_e(2) > 1 && x_e(1) < m && x_e(2) < n)
            I2(i,j,:) = I(x_e(1),x_e(2),:);
        end
    end
end


end

