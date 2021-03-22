function N = homography_fit(P, H, tau)
%% pass points through homography
p1 = [P(:,1:2), ones(size(P,1), 1)].';
p2hat = (H * p1).';
p2hat = p2hat(:,1:2);

%% compare to found correspondences
p2 = P(:,3:4);
err = sum((p2 - p2hat).^2, 2);
N = err < tau;
end