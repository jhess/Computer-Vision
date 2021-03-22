function [u, v] = opticalFlowLK(I1, I2, W, t)
warning('off','all')
%OPTICALFLOWJK Summary
% Usage: corners = harris(I, sigma, threshold)
% Inputs:
%   I1  [m x n]  Input image 1 to be processed
%   I2  [m x n]  Input image 2 to be processed
%   W    Window size specifying area of neighborhood around each pixel
%   t   threshold for comparing to eigenvalues of left matrix
% Outputs:
%   u, v The output velocity vector 

% Convert to grayscale
if size(I1, 3) == 3
    I1 = rgb2gray(I1);
end
I1 = mat2gray(I1);
if size(I2, 3) == 3
    I2 = rgb2gray(I2);
end
I2 = mat2gray(I2);

% Convert to double floating precision
I1 = double(I1);
I2 = double(I2);

%gaussian filter smoothing
gaussian_filter = fspecial('gaussian', [3,3], 1);
I1_smoothed = imfilter(I1, gaussian_filter);
I2_smoothed = imfilter(I2, gaussian_filter);

%compute gradients at pixels
hx = 1/12.*[-1,8,0,-8,1];
%hx = [1 0 -1];
hy = hx';

% compute the gradient/derivative with respect to x of smoothed I1
Ix = conv2(I1_smoothed, hx, 'same');

% compute the gradient with respect to y of smoothed I1
Iy = conv2(I1_smoothed, hy, 'same');

%compute the temporal gradient
It = I2_smoothed - I1_smoothed;

% Compute quadratic terms of the gradients/derivatives
%Ix2 = Ix.^2;  
%Iy2 = Iy.^2;
%Ixy = Ix.*Iy;
Ix2 = conv2(Ix.^2, ones(W),'same');
Iy2 = conv2(Iy.^2, ones(W),'same');
Ixy = conv2(Ix.*Iy, ones(W),'same');

Ixt = conv2(Ix.*It, ones(W),'same');
Iyt = conv2(Iy.*It, ones(W),'same');

%get size of image
s = size(I1);

u = zeros(s);
v = zeros(s);

%use input window, W
hs = floor(W/2);
for i = 1:s(1)
    for j =1:s(2)
        left = j - hs; 
        right = j + hs;
        top = i - hs; 
        bottom = i + hs;
        if(left<=0)
            left=1; 
        end
        if(right>s(2))
            right=s(2); 
        end
        if(top<=0)
            top = 1; 
        end
        if(bottom>s(1))
            bottom=s(1); 
        end
        %region
        region = (right - left + 1)*(bottom - top + 1);
        %region = I1(top:bottom,left:right);
        %ix = Ix(top:bottom,left:right);
        %iy = Iy(top:bottom,left:right);
        
        % left matrix, corner response matrix of quadratic gradients
        A = [Ix2(i,j) Ixy(i,j); Ixy(i,j) Iy2(i,j)]/region;

        %right matrix, temporal quadratics
        B = (-1).*[Ixt(i,j); Iyt(i,j)]/region;
        
        %compute the rank of the gradient left matrix
        r = rank(A);
        %compute the eigenvalues of corner reponse left matrix
        lambda = eig(A);
        %if min eigenvalue is greater than threshold
        if (min(lambda)>t)
            %hitMap(i,j) = 1;
            %solve system of linear equations for velocity vector
            UV = A\B;
            u(i,j) = UV(1);
            v(i,j) = UV(2);
        end
        
    end
end
        

end

