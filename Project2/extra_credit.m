%%
% Extra Credit
insert = imread('belichick.jpg');
hallway = imread('DanaHallWay2/DSC_0285.JPG');
figure(8)
imshow(hallway)
hold on
% input the 4 corner points where to insert, clockwise, starting from top left
[x,y] = ginput(4); 
plot(x, y, 'gx')
hold off
% matrix map of background hallway image
hallway_m = [1 1;size(insert,2) 1;size(insert,2) size(insert,1);1 size(insert,1)];
% create a matrix map of the input corners
insert_m = [x y];
% find the homography matrix
[h, ~, ~] = estimateGeometricTransform(hallway_m, insert_m, 'projective');
warped_mosaic = imwarp(insert,h);
disp(insert_m);
% replace the ginput surface with insert image, and shift the ginput limits
insert_bounds = [insert_m(:,1)-min(insert_m(:,1))+1 insert_m(:,2)-min(insert_m(:,2))+1];
% create preliminary mask of zeroes for warped image
mask = zeros(ceil(fliplr(max(insert_bounds))));
for i = 1:size(mask,1)
    copy_warped = repmat(i,1,size(mask,2));
    size_warped = 1:size(mask,2);
    % check if warped image are within the bounds of insert image bounds first
    % two columns, x and y
    % populate each row of the mask
    mask(i,:) = inpolygon(copy_warped, size_warped, insert_bounds(:,2), insert_bounds(:,1));
end
% overlay image
coord = round(min(insert_m));%top left
for i = 1:size(mask,1)
    for j = 1:size(mask,2)
        if mask(i,j)==1
            % set the coordinates of hallway image within our bounds to the
            % inserted input image to overlay onto
            hallway(coord(2) + i-1, coord(1) + j-1, :) = warped_mosaic(i, j, :);
        end
    end
end
%show the final background hallway image with the insert picture overlayed
%on it
imshow(hallway)
hold off