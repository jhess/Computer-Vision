function H = homography_ransac(P, tau, max_iter, alpha)
%% begin the RANSAC iterations
best_match = 0;
good_enough = @(N) N > alpha * size(P,1);
for n = 1:max_iter
    %% pick points
    idx = randperm(size(P,1), 4);
    Phat = P(idx, :);
    %% find a homography
    Hhat = homography(Phat);
    
    %% check the number of matching points
    fit = homography_fit(P, Hhat, tau);
    Nhat = sum(fit);
    
    %% check for termination conditions
    if good_enough(Nhat)
        H = homography(P(fit, :));
        return
    end
    
    %% update best match
    if Nhat > best_match
        best_match = Nhat;
        best_idx = fit;
    end
end
if best_match == 0
    error('No match found. tau or alpha may be too restrictive.');
end
H = homography(P(best_idx,:));
end