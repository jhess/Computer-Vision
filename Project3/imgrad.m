function [Ix, Iy] = imgrad(I, hx, hy)
% filters
if nargin == 1
    hx = [-1, 0, 1;...
          -2, 0, 2;...
          -1, 0, 1];
end
if nargin < 3
    hy = hx.';
end

% filter the image
Ix = imfilter(I, hx);
Iy = imfilter(I, hy);
end