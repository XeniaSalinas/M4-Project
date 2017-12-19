function cropped_I = crop_transformed_image(I)
% CROP_TRANSFORMED_IMAGE Crops an image transformed through successive
% projective transformations (with apply_H)

% Ensure image has datatype uint8
I_uint = uint8(I);

% Compute mask and bounding box
mask = I_uint(:,:,1) > 0 | I_uint(:,:,2) > 0 | I_uint(:,:,3) > 0;
im_props = regionprops(mask, 'BoundingBox');
bb = round(im_props.BoundingBox);
x = bb(1);
y = bb(2);
w = bb(3);
h = bb(4);

cropped_I = I(y:y+h-1,x:x+w-1,:);

end

