function [corners] = harris(I, sigma, threshold)
warning('off','all')
% Usage: corners = harris(I, sigma, threshold)
% Inputs:
%   I  [m x n]  Input image to be processed
%   sigma    Standard deviation of the Gaussian filter
%              used for smoothing the derivatives
%   threshold  Corner response function threshold. Values of the corner
%              response function below this threshold are neglected
% Outputs:
%   corners  [n x 3]  Each row of 3-column matrix represents
%     one corner point. The first two columns describe the position
%     of the corner in im_inp (first column = y , second column = x). 
%     The third column contains the value of the corner response function.

% Conversion of the limiting variables to doubles
sigma = double(sigma);

% Convert im to double precision
I = double(I);

image(I);

% Firstly, the image is smoothed with a simple gaussian filter

sig = 1;
window1 = floor(sig*3);
x = -window1 : 1 : window1; 
gaussian1DFilter = exp(-(x).^2 / (2*sig^2)) / (sig*sqrt(2*pi));
   
% smoothed image
I_smoothed = conv2(gaussian1DFilter, gaussian1DFilter', I, 'same' );

% Compute image derivatives by convolution with kernel
% [1 0 -1] in both x and y directions

hx = [1 0 -1];
hy = hx';
%hx = [1, 0, -1; 1, 0, -1; 1, 0, -1];
%hy = [1, 0, -1; 1, 0, -1; 1, 0, -1].';

% compute the gradient/derivative with respect to x 
Ix = conv2(I_smoothed, hx, 'same');

% compute the gradient with respect to x at the border pixels 
Ix(:,1)   = 2*(I_smoothed(:,2)  -I_smoothed(:,1));
Ix(:,end) = 2*(I_smoothed(:,end)-I_smoothed(:,end-1));

% compute the gradient with respect to y 
Iy = conv2(I_smoothed, hy, 'same');

% compute the gradient with respect to y at the border pixels 
Iy(1,:)   = 2*(I_smoothed(2,:)  -I_smoothed(1,:));
Iy(end,:) = 2*(I_smoothed(end,:)-I_smoothed(end-1,:));

% Precompute quadratic terms of the derivatives.
Ix2 = Ix.^2;  
Iy2 = Iy.^2;
Ixy = Ix.*Iy;

% We then smooth the derivatives/gradients with another gaussian filter
% with a standard deviation of the second parameter input, sigma

window2 = floor(3*sigma);
x = -window2:1:window2;
gaussian_gradient_filter = exp(-(x).^2/(2*sigma^2)) / (sigma*sqrt(2*pi));

Ix2 = conv2( gaussian_gradient_filter, gaussian_gradient_filter', Ix2, 'same');
Iy2 = conv2( gaussian_gradient_filter, gaussian_gradient_filter', Iy2, 'same');
Ixy = conv2( gaussian_gradient_filter, gaussian_gradient_filter', Ixy, 'same');

% Set sensitivity factor of the response function,
% used to detect sharp corners (usually 0.04 <= k <= 0.06)

k = 0.04;

% Corner response, r
% calculate the determinant and trace of C matrix elements
% and then calculate r

detC = (Ix2.*Iy2 - Ixy.^2);
traceC = Ix2 + Iy2;
r = detC - k * traceC.^2;

figure(3)
imagesc(r), colormap(hot(256)), axis image, colorbar,
title('Corner Response values')

r_thresholded = r .* (r > threshold);

figure(4)
imagesc(r_thresholded.^(1/2)), colormap(hot(256)), axis image, colorbar,
title('Thresholded Corner Response Output');

%  Non-maximal suppression
%  loop over corner points and only keep those where the corner response 
%  is a peak/local maximum in a 3 x 3 window. All other pixels
%  are set to zero if they do not pass the threshold.

r_final = r_thresholded;

% create a temporary pad and threshold
pad = zeros( size(r_thresholded) + [2 2] );
pad( 2: end-1, 2 : end-1) = r_thresholded;

r_thresholded = pad;

% loop over non-zero pixels that passed the threshold
% pixels that don't pass threshold are set to 0
[idx, idy] = find(r_thresholded>0);

for k = 1:length(idx)
  i_x = idx(k);  
  i_y = idy(k);
  % compare to maximum of the 8-pixel neighborhood
  if r_thresholded(i_x, i_y) ~= max(max(r_thresholded(i_x-1 : i_x+1,...
          i_y - 1 : i_y+1))) 
  % use -1 here because of the padding
  r_final(i_x-1,i_y-1) = 0;   
  end
end

%  Output 

%  corner points are represented by a matrix with the 
%  coordinate values x and y as the first two and the third
%  is the value of the corner response 

[cornersx, cornersy, response] = find(r_final);

corners = [cornersx cornersy response];

return; 

end

