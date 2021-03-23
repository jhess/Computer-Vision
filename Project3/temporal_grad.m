function [It] = temporal_grad(I)
% filter impulse response
h = [-1, 1];
h = permute(h, [1, 3, 2]);

% filter the video
It = imfilter(I, h, 'full');
It = It(:,:,1:size(I,3)); % dependent on filter length
end