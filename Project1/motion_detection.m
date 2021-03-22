clear all
close all
clc
imtool close all

image = 'EnterExitCrossingPaths2cor%04d.jpg';
im_start = 0;
im_end = 484;

%Convert images into grayscale and read them as a movie
num_images = im_end-im_start+1;
for i = im_start:im_end
    images(:,:,i+1) = rgb2gray(imread(sprintf(image,i)));
end

%Create the temporal filter
filter_matrix = [-1 0 1];
f_coeff = 0.5;
temp_filter = f_coeff*filter_matrix;

%Get depth of temporal filter
depth = (length(temp_filter) - 1)/2;

%initialize the filtered image with zeros
filtered_image = zeros(size(images(:,:,1),1),size(images(:,:,1),2),num_images);
%create the mask
image_mask = filtered_image;


for i = 1+depth:num_images-depth
    %initialize the derivative pixel intensity values
    derivative_values = zeros(size(images(:,:,1)));
    for l = 1:length(temp_filter)
        %compute the derivatives
        derivative_values = derivative_values + double(images(:,:,i-depth+l-1))*temp_filter(l);
    end
    %get absolute value of the derivatives of the image values
    filtered_image(:,:,i) = abs(derivative_values);
    %imshow(filtered_image(:,:,i));
end

%set the threshold of temporal derivative intensity values
threshold = 18;
%apply the thresholding to the mask
image_mask(filtered_image >= threshold) = 1;
image_mask(filtered_image < threshold) = 0;

motion_result = images;
%set the result of mask==1 to maximum grayscale pixel value
motion_result(image_mask==1) = 255;

%create the video of motion tracking
implay(motion_result)

%%
%3x3 box filter
b_size = 3;
for i = 1:num_images
    images_box(:,:,i) = imfilter(images(:,:,i), fspecial('average',[b_size,b_size]));
end

temp_filter = 0.5*[-1 0 1];
filter_depth = (length(temp_filter)-1)/2;
filtered_image = zeros(size(images_box(:,:,1),1), size(images_box(:,:,1),2),num_images);
image_mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    derivative_values = zeros(size(images_box(:,:,1)));
    for l = 1:length(temp_filter)
        derivative_values = derivative_values + double(images_box(:,:,i-filter_depth+l-1))*temp_filter(l);
    end
    filtered_image(:,:,i) = abs(derivative_values);
end

%set the threshold of temporal derivative intensity values
threshold = 18;
%apply the thresholding to the mask
image_mask(filtered_image >= threshold) = 1;
image_mask(filtered_image < threshold) = 0;

motion_result = images_box;
%set the result of mask==1 to maximum grayscale pixel value
motion_result(image_mask==1) = 255;

%create the video of motion tracking
implay(motion_result)

%%
%5x5 box filter
b_size = 5;
for i = 1:num_images
    images_box(:,:,i) = imfilter(images(:,:,i),fspecial('average',[b_size,b_size]));
end

temp_filter = 0.5*[-1 0 1];
filter_depth = (length(temp_filter)-1)/2;
filtered_image = zeros(size(images_box(:,:,1),1),size(images_box(:,:,1),2),num_images);
image_mask = filtered_image;

for i = 1+filter_depth:num_images-filter_depth
    derivative_values = zeros(size(images_box(:,:,1)));
    for l = 1:length(temp_filter)
        derivative_values = derivative_values + double(images_box(:,:,i-filter_depth+l-1))*temp_filter(l);
    end
    filtered_image(:,:,i) = abs(derivative_values);
end
%set the threshold of temporal derivative intensity values
threshold = 18;
%apply the thresholding to the mask
image_mask(filtered_image >= threshold) = 1;
image_mask(filtered_image < threshold) = 0;

motion_result = images_box;
%set the result of mask==1 to maximum grayscale pixel value
motion_result(image_mask==1) = 255;

%create the video of motion tracking
implay(motion_result)