clear all
close all
clc
imtool close all

%I1 = imread('corridor/bt.000.png');
%I2 = imread('corridor/bt.001.png');

I = zeros(100, 100, 2);
I(40:60, 40:60, 1) = 1;
I(:, :, 2) = circshift(I(:,:,1), [0, 10]);

I1 = I(:, :, 1);
I2 = I(:, :, 2);

% Convert to grayscale
if size(I1, 3) == 3
    I1 = rgb2gray(I1);
end
I1 = mat2gray(I1);
if size(I2, 3) == 3
    I2 = rgb2gray(I2);
end
I2 = mat2gray(I2);

% Segment into Two-level Pyramids
i1 = cell(3,1);
i1{3} = I1;
i1{2} = impyramid(I1, 'reduce');
i1{1} = impyramid(i1{2}, 'reduce');


i2 = cell(3,1);
i2{3} = I2;
i2{2} = impyramid(I2, 'reduce');
i2{1} = impyramid(i2{2}, 'reduce');

t = 0.2;
windowSize = [15,30,100];

% optical flow
for w = 1:3
    wsize = windowSize(w);
    [u, v] = opticalFlowLK(I1, I2, wsize, t);

    [x, y] = meshgrid(1:10:size(I1,1), size(I1,1):-10:1);
    
    qu = u(1:10:size(I1,1), 1:10:size(I1,2));
    qv = v(1:10:size(I1,1), 1:10:size(I1,2));
    
    subplot(2,3,w), quiver(x,y, qu, -qv,'linewidth', 1),axis([1,size(I1,1),1,size(I1,2)]), title(sprintf('windowSize=%d',wsize));
end