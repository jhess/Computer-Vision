clc; close all; clear;

% parameters
tau = 100;
max_iter = 100;
alpha = 0.9;

% read in 2 images
hallway_original1 = imread('DanaHallWay2/DSC_0286.JPG'); 
hallway_original2 = imread('DanaHallWay2/DSC_0287.JPG');

I1 = rgb2gray(hallway_original1);
I2 = rgb2gray(hallway_original2);

% apply Harris corner detection
c1 = harris(I1, 4, 25000);
c2 = harris(I2, 4, 25000);

% show detected corners
figure(5)
imshow(hallway_original1)
hold on
plot(c1(:,2),c1(:,1),'r*')
hold off
figure(6)
imshow(hallway_original2)
hold on
plot(c2(:,2),c2(:,1),'r*')
hold off

% Normalized Cross Correlation

% corrs are the correlation points
corrs = [];
for i = 1:size(c1, 1)
    
    y1 = c1(i, 1);
    x1 = c1(i, 2);
    % p_size = patch size, patch length = 2*patch_size, radius of patch
    p_size = 10;
    
    % ignore border corners
    if x1 <= p_size || y1 <= p_size || x1 >= size(I1,2)-p_size || y1 >= size(I1,1)-p_size
        continue
    end
    
    %creates patch around each corner for both images
    %p1 is patch 1
    p1 = I1(y1-p_size : y1+p_size,x1-p_size : x1+p_size);
    
    % initialize NCC to a matrix of zeroes
    ncc = zeros(1,size(c2,1));
    
    for j = 1:size(c2,1)
        y2 = c2(j, 1);
        x2 = c2(j, 2);

        if x2 <= p_size || y2 <= p_size || x2 >= size(I2,2)-p_size || y2 >= size(I2,1)-p_size
            continue
        end
        
        % p2 is patch 2
        p2 = I2(y2-p_size : y2+p_size,x2-p_size : x2+p_size);
        
        % correlation coefficients in a matrix obtained by 
        % 2D normalized cross correlation**
        corr_coeffs = normxcorr2(p1, p2);
        
        % obtain center value of ncc matrix for correlation value between
        % 2 patches
        ncc(j) = corr_coeffs(1+2*p_size, 1+2*p_size);
        
        % find the peak in the cross correlation
        if ncc(j) == max(ncc)
            c_index = j;
        end
    end
    
    % thresholding
    ncc_threshold = 0.90;
    if ncc(c_index) > ncc_threshold
        % set the correlation points after if the peak cross correlation
        % passes the threshold
        corrs = [corrs;x1 y1 c2(c_index,2) c2(c_index,1)];
    end
end

figure(7)
imshow([hallway_original1;hallway_original2])
hold all
plot(c1(:,2), c1(:,1), 'rx');
plot(c2(:,2), c2(:,1) + size(I1,1), 'bx');
for i = 1:size(corrs,1)
    plot([corrs(i,1), corrs(i,3)], ...
        [corrs(i,2), corrs(i,4)+size(I1,1)], 'g');
end
hold off

% find a homography
H = homography_ransac(corrs, tau, max_iter, alpha);
disp('Homography:');
disp(H);

%H2 = [H [0; 0; 1]];
%disp(H2);
cf1 = corrs(:,1:2);
cf2 = corrs(:,3:4);
% warp the image
T = affine2d(H.');
%T = estimateGeometricTransform(cf1,cf2,'projective');
I2 = imwarp(I1, T);

% mosaic
Iout = warp_img(I2, H, I1);

% plot
figure;
imagesc(I1);
figure;
imagesc(I2);
figure;
imagesc(Iout);
% figure;
% imagesc(R1);
% figure;
% imagesc(R2);

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
            hallway(coord(2)+i-1,coord(1)+j-1,:) = warped_mosaic(i,j,:);
        end
    end
end
imshow(hallway)
hold off