function H = homography(P)
%% normalize
P1 = [P(:, 1:2), ones(size(P,1), 1)];
H1 = normalization_transform(P1);
P1 = (H1 * P1.').';
P2 = [P(:, 3:4), ones(size(P,1), 1)];
H2 = normalization_transform(P2);
P2 = (H2 * P2.').';
P = [P1(:,1:2), P2(:,1:2)];

%% form the system of linear equations
A = zeros(2 * size(P, 1), 8);
b = zeros(2 * size(P, 1), 1);
for n = 1:size(P, 1)
    %% extract the points
    x1 = P(n, 1);
    y1 = P(n, 2);
    x2 = P(n, 3);
    y2 = P(n, 4);
    %% A matrix entries
    A(2*n, :)     = [x1, y1, 1, 0, 0, 0, -x1*x2, -y1*x2];
    A(2*n + 1, :) = [0, 0, 0, x1, y1, 1, -x1*y2, -y1*y2];
    %% b matrix entries
    b(2*n, 1)     = x2;
    b(2*n + 1, 1) = y2;
end

%% solve
h = pinv(A) * b;

%% determine the transformation matrix
h = [h; 1];
H = reshape(h, 3, 3).';

%% unnormalize
H = H2^-1 * H * H1;
end