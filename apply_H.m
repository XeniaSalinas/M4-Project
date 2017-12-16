function I2 = apply_H(I, H)
%APPLY_H Transform the input image with the given homography

%Compute the size of the transformed image
[n,m,c] = size(I);
lu_p = [1;1;1];
ru_p = [1;m;1];
ld_p = [n;1;1];
rd_p = [n;m;1];

lu2_p = H * lu_p;
ru2_p= H * ru_p;
ld2_p = H * ld_p;
rd2_p = H * rd_p;

lu2_e = [lu2_p(1)/lu2_p(3), lu2_p(2)/lu2_p(3)];
ru2_e = [ru2_p(1)/ru2_p(3), ru2_p(2)/ru2_p(3)];
ld2_e = [ld2_p(1)/ld2_p(3), ld2_p(2)/ld2_p(3)];
rd2_e = [rd2_p(1)/rd2_p(3), rd2_p(2)/rd2_p(3)];

n2 = ceil(abs(max([lu2_e(1), ru2_e(1), ld2_e(1), rd2_e(1)]) - min([lu2_e(1), ru2_e(1), ld2_e(1), rd2_e(1)]))) + 1;
m2 = ceil(abs(max([lu2_e(2), ru2_e(2), ld2_e(2), rd2_e(2)]) - min([lu2_e(2), ru2_e(2), ld2_e(2), rd2_e(2)]))) + 1;

%Compute the position rectification vector
origin_n = min([lu2_e(1), ru2_e(1), ld2_e(1), rd2_e(1)]);
origin_m = min([lu2_e(2), ru2_e(2), ld2_e(2), rd2_e(2)]);

%Compute the homography transfromation
I2 = zeros(n2,m2,c);

for i = 1:n2
    for j = 1:m2
        x2_e = [i; j] + [origin_n; origin_m];
        x2_p = H\[x2_e;1]; %A\b for inv(a)*b
        x_e = round([x2_p(1)/x2_p(3), x2_p(2)/x2_p(3)]); %TODO: Use interp2
        if (x_e(1) > 1 && x_e(2) > 1 && x_e(1) < n && x_e(2) < m)
            I2(i,j,:) = I(x_e(1),x_e(2),:);
        end
    end
end


end

