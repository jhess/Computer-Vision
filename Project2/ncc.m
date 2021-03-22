%% Normalized Cross Correlation

% corrs are the correlation points
corrs = [];
for i = 1:size(c1, 1)
    
    y1 = c1(i, 1);
    x1 = c1(i, 2);
    % p_size = patch size, patch length = 2*patch_size, radius of patch
    p_size = 10;
    
    % ignore border corners
    if x1 <= p_size || y1 <= p_size || x1 >= size(I1,2)-p_size || y1 >= size(I1,1)-p_size
        continue
    end
    
    %creates patch around each corner for both images
    %p1 is patch 1
    p1 = I1(y1-p_size : y1+p_size,x1-p_size : x1+p_size);
    
    % initialize NCC to a matrix of zeroes
    ncc = zeros(1,size(c2,1));
    
    for j = 1:size(c2,1)
        y2 = c2(j, 1);
        x2 = c2(j, 2);

        if x2 <= p_size || y2 <= p_size || x2 >= size(I2,2)-p_size || y2 >= size(I2,1)-p_size
            continue
        end
        
        % p2 is patch 2
        p2 = I2(y2-p_size : y2+p_size,x2-p_size : x2+p_size);
        
        % correlation coefficients in a matrix obtained by 
        % 2D normalized cross correlation**
        corr_coeffs = normxcorr2(p1, p2);
        
        % obtain center value of ncc matrix for correlation value between
        % 2 patches
        ncc(j) = corr_coeffs(1+2*p_size, 1+2*p_size);
        
        % find the peak in the cross correlation
        if ncc(j) == max(ncc)
            c_index = j;
        end
    end
    
    % thresholding
    ncc_threshold = 0.90;
    if ncc(c_index) > ncc_threshold
        % set the correlation points after if the peak cross correlation
        % passes the threshold
        corrs = [corrs;x1 y1 c2(c_index,2) c2(c_index,1)];
    end
end

figure(7)
imshow([hallway_original1;hallway_original2])
hold all
plot(c1(:,2), c1(:,1), 'rx');
plot(c2(:,2), c2(:,1) + size(I1,1), 'bx');
for i = 1:size(corrs,1)
    plot([corrs(i,1), corrs(i,3)], ...
        [corrs(i,2), corrs(i,4)+size(I1,1)], 'g');
end
hold off