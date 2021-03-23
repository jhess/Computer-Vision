function Iout = warp_img(I, H, Iref)
% handle exceptional cases for the reference image
if nargin == 2
    Iref = [];
end
if numel(Iref) == 0
    Iref = zeros(2);
end

% get points for the images
[y, x] = ndgrid(0:size(I,1)-1, 0:size(I,2)-1);
P = [x(:), y(:), ones(numel(x), 1)];
[yref, xref] = ndgrid(0:size(Iref,1)-1, 0:size(Iref,2)-1);
Pref = [xref(:), yref(:), ones(numel(xref), 1)];

% transform points from the image
Hinv = H ^ -1;
Phat = apply_homography(P, Hinv);

% find the size for the output image
min_vals = min(Phat(:,1:2),[],1);
min_vals = min([Pref(:,1:2); min_vals], [], 1);
min_vals = floor(min_vals);
max_vals = max(Phat(:,1:2), [], 1);
max_vals = max([Pref(:,1:2); max_vals], [], 1);
max_vals = ceil(max_vals);
out_sz = max_vals - min_vals;
out_sz = [out_sz(2), out_sz(1)];

% define the output image and coordinate system
Iout = zeros(out_sz);
count = zeros(out_sz);
[yout, xout] = ndgrid(0:size(Iout,1)-1, 0:size(Iout,2)-1);
yout = yout + min_vals(2);
xout = xout + min_vals(1);

% warp the reference image
Ioutref = interp2(Iref, xout, yout);
region = not(isnan(Ioutref));
count(region) = count(region) + 1;
Iout(region) = Iout(region) + Ioutref(region);

% reverse transform the output image coordinates
Pout = [xout(:), yout(:), ones(numel(xout),1)];
Pout = apply_homography(Pout, H);
xout = reshape(Pout(:,1), size(xout));
yout = reshape(Pout(:,2), size(yout));

% warp the main image
Ioutmain = interp2(I, xout, yout);
region = not(isnan(Ioutmain));
count(region) = count(region) + 1;
Iout(region) = Iout(region) + Ioutmain(region);

% average
region = count ~= 0;
Iout(region) = Iout(region) ./ count(region);
% figure;
% imagesc((0:size(Iout,2)) + min_vals(2), (0:size(Iout,1)) + min_vals(1), Iout);
% colormap('gray');

end