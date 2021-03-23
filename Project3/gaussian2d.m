function h = gaussian2d(M, N, cov)
% get the spatial coordinates
x1 = (0.5:M-0.5) - M/2;
x2 = (0.5:N-0.5) - N/2;
[X1, X2] = ndgrid(x1, x2);
X = [reshape(X1, [], 1), reshape(X2, [], 1)];

% evaluate the Gaussian
h = exp(-0.5 * dot(X * cov^-1, X, 2)) / sqrt(det(2 * pi * cov));
h = reshape(h, M, N);

% make sum to one
h = h / sum(h(:));
end