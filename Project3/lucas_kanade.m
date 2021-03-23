function [u, v] = lucas_kanade(I, w, tau)
if nargin < 3
    tau = 1e-4;
end
N = 5;
I = imfilter(I, box2d(N, N));
% compute spatial gradient
hx = -[-1, 8, 0, -8, 1] / 12;
[Ix, Iy] = imgrad(I, hx);

% compute temporal gradient
It = temporal_grad(I);

% smooth the gradients
% h = [1, 4, 6, 4, 1] / 16;
% h = ones(1, 7) / 7;
% Ix = imfilter(imfilter(Ix, h), h.');
% Iy = imfilter(imfilter(Iy, h), h.');
% It = imfilter(imfilter(It, h), h.');

% compute gradient products
Ix2 = Ix.^2;
Iy2 = Iy.^2;
Ixy = Ix .* Iy;
Ixt = Ix .* It;
Iyt = Iy .* It;


% apply neighborhood
Ix2 = imfilter(Ix2, w, 'same');
Iy2 = imfilter(Iy2, w, 'same');
Ixy = imfilter(Ixy, w, 'same');
Ixt = imfilter(Ixt, w, 'same');
Iyt = imfilter(Iyt, w, 'same');

% solve system of equations
u = zeros(size(I));
v = zeros(size(I));
e1 = zeros(size(I));
e2 = zeros(size(I));
for m = 1:size(I, 1)
    for n = 1:size(I,2)
        for p = 1:size(I,3)
            A = [Ix2(m,n,p), Ixy(m,n,p); Ixy(m,n,p), Iy2(m,n,p)];
            b = [-Ixt(m,n); -Iyt(m,n)];
            e = eig(A);
            e1(m,n,p) = e(1);
            e2(m,n,p) = e(2);
            if all(e > tau)
                % matrix is non-singular. just solve the system of
                % equations
%                 x = A^-1 * b;
%                 u(m,n,p) = x(2);
%                 v(m,n,p) = x(1);
                detA = det(A);
                v(m,n,p) = det([b, A(:,2)]) / detA;
                u(m,n,p) = det([A(:,1), b]) / detA;
            elseif any(e > tau)
                % matrix is singular but we are still on an edge. estimate
                % the optical flow as the normal velocity vector
                grad = [Ix(m,n,p), Iy(m,n,p)];
                grad_norm = norm(grad);
                s = -It(m, n, p) / grad_norm;
                if grad_norm < 1e-3
                    norm_velocity = [0, 0];
                else
                    norm_velocity = s * grad / grad_norm;
                end
                u(m,n,p) = norm_velocity(1);
                v(m,n,p) = norm_velocity(2);
            end
        end
    end
end

% figure;
% imdisp(e1(:,:,1));
% figure;
% imdisp(e2(:,:,1));
% figure;
% imdisp(e1(:,:,2));
% figure;
% imdisp(e2(:,:,2));
% figure;
% imdisp(It(:,:,1));
figure;
imdisp(It(:,:,2));
% figure;
% imdisp(Ix(:,:,1));
figure;
imdisp(Ix(:,:,2));
% figure;
% imdisp(Iy(:,:,1));
figure;
imdisp(Iy(:,:,2));
figure;
imdisp(Ix2(:,:,2));
figure;
imdisp(Iy2(:,:,2));
figure;
imdisp(Ixy(:,:,2));
figure;
imdisp(Ixt(:,:,2));
figure;
imdisp(Iyt(:,:,2));
end