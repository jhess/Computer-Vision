%% load video/extract frames
clear all
close all

video_path = 'C:\Users\justhess\Documents\MATLAB\Project4\data\tiger1\imgs';

contents = dir(video_path);
names = {};
for k = 1:numel(contents)
    name = contents(k).name;
    if ~strcmp(name, '.') && ~strcmp(name, '..')
        names{end+1} = name;
    end
end
num_frames = numel(names);

frames = [];
for i = 1:num_frames
    frames(:,:,i) = rgb2gray(imread(fullfile(video_path,names{i})));
end
frames = uint8(frames); %rgb friendly

%% setup

% parameters according to the paper
padding = 1;					%extra area surrounding the target
output_sigma_factor = 1/16;		%spatial bandwidth (proportional to target)
sigma = 0.2;					%gaussian kernel bandwidth
lambda = 1e-2;					%regularization
interp_factor = 0.075;			%linear interpolation factor for adaptation

%pos = box_position(1:2)+floor(box_position(3:4)/2);
start_pos = [65 135];%fliplr(pos); %hardcode positions for consistent presentation
%target_sz = box_position(3:4);
target_sz = [42 38];%fliplr(target_sz); %hardcode positions for consistent presentation

%window size, taking padding into account
sz = floor(target_sz * (1 + padding));

%desired output (gaussian shaped), bandwidth proportional to target size
output_sigma = sqrt(prod(target_sz)) * output_sigma_factor;
[rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
y = exp(-0.5 / output_sigma^2 * (rs.^2 + cs.^2));
yf = fft2(y);

%store pre-computed cosine window
cos_window = hann(sz(1)) * hann(sz(2))';


%% iterate over frames

use_kalman = true; %t for kalman based occlusion recovery, f for only tracking

pos = start_pos;

time = 0;  %to calculate FPS
positions = zeros(num_frames, 2);  %to calculate precision, these are centers
rng(1)
PSR_values = [];
pos_values = [];
if use_kalman
    occlusion_frames = [];
end

for i = 1:num_frames
    
    %for screenshot purposes
    occlude_flag = false;
    
    %load current image from frame stack
    img = frames(:,:,i);
    
    %extract and pre-process subwindow of full image
    x = get_subwindow(img, pos, sz, cos_window);
    
    %default box color
    rect_color = 'g';
    
    %do processing for all frames after first
    if i > 1
        %calculate response of the classifier at all locations
        k = dense_gauss_kernel(sigma, x, z);
        response = real(ifft2(alphaf .* fft2(k)));   %(Eq. 9)
        
        %target location is at the maximum response
        [row, col] = find(response == max(response(:)), 1);
        pos = pos - floor(sz/2) + [row, col];
        pos_values(i,:) = pos;
        
        if use_kalman
            %calculate PSR (peak sidelobe ratio)
            gmax = response(row,col);
            
            window_size = round(min(sz)/7);
            %window_size = 11;
            sq_radius = (window_size-1)/2;
            %sidelobe area is target area - peak window area
            sidelobe = zeros(1,numel(response)-window_size^2);
            % we can store this as 1d since spatial properties dont matter
            % here
            idx = 1;
            for r = 1:size(response,1)
                for c = 1:size(response,2)
                    %sidelobe consists of values outside target window around max response
                    if (r < row - sq_radius || r > row + sq_radius || c < col - sq_radius || c > col + sq_radius)
                        sidelobe(idx) = response(r,c);
                        idx = idx+1;
                    end
                end
            end
            sidelobe_mean = mean(sidelobe);
            sidelobe_std = std(sidelobe);

            PSR = (response(row,col) - sidelobe_mean)/sidelobe_std;
            PSR_values(i) = PSR;
            
            
            %check for occlusion based on PSR
            %PSR_threshold = prctile(PSR_values,1);
            PSR_threshold = 4.6;
            if PSR < PSR_threshold
                occlusion_frames = [occlusion_frames i];
                occlude_flag = true;
                rect_color = 'r';
                %pause(0.5)
                fprintf('Occlusion detected in frame %d, PSR = %g\n',i,PSR)
                
                %kalman filter prediction
                Nsamples = 10;
                dt = 1; %1 frame
                t=0:dt:dt*Nsamples;
                Vtrue_x = (pos_values(end-1,2)-pos_values(end-Nsamples,2))/(Nsamples-1);
                Vtrue_y = (pos_values(end-1,1)-pos_values(end-Nsamples,1))/(Nsamples-1);
                Xinitial = pos_values(end-Nsamples,2);
                Yinitial = pos_values(end-Nsamples,1);
                Xtrue = Xinitial + Vtrue_x * t;
                Ytrue = Yinitial + Vtrue_y * t;
                Xk_prev = [0;.5*Vtrue_x];
                Yk_prev = [0;.5*Vtrue_y];
                Xk = [];
                Yk = [];
                Phi = [1 dt;0  1];
                sigma_model = 1;
                % P determines to weight measurement or estimate more
                P = [sigma_model^2 0;0 sigma_model^2];
                % Q = process noise covar
                Q = [0 0; 0 0];
                % M is measurement [X V], v is 0 since we only measure x
                M = [1 0];
                sigma_meas = 1;
                R = sigma_meas^2;
                
                Xk_buffer = zeros(2,Nsamples+1);
                Yk_buffer = zeros(2,Nsamples+1);
                Xk_buffer(:,1) = Xk_prev;
                Yk_buffer(:,1) = Yk_prev;
                Z_buffer_x = zeros(1,Nsamples+1);
                Z_buffer_y = zeros(1,Nsamples+1);
                for h=1:Nsamples
                    
                    % Z is the measurement vector. In our
                    % case, Z = TrueData + uncertainty from occlusion
                    % (modeled as gaussian)
                    Z_x = Xtrue(h+1)+sigma_meas*randn;
                    Z_y = Ytrue(h+1)+sigma_meas*randn;
                    Z_buffer_x(h+1) = Z_x;
                    Z_buffer_y(h+1) = Z_y;
                    
                    % Kalman iteration
                    P1 = Phi*P*Phi' + Q;
                    S = M*P1*M' + R;
                    
                    % K is Kalman gain. If K is large, more weight goes to the measurement.
                    % If K is low, more weight goes to the model prediction.
                    K = P1*M'*inv(S);
                    P = P1 - K*M*P1;
                    
                    Xk = Phi*Xk_prev + K*(Z_x-M*Phi*Xk_prev);
                    Yk = Phi*Yk_prev + K*(Z_y-M*Phi*Yk_prev);
                    Xk_buffer(:,h+1) = Xk;
                    Yk_buffer(:,h+1) = Yk;
                    
                    % For the next iteration
                    Xk_prev = Xk;
                    Yk_prev = Yk;
                end
                pos = [Yk_buffer(1,end) Xk_buffer(1,end)];
                %update position with kalman filter estimate
                pos_values(i,2) = round(pos(2));
                pos_values(i,1) = round(pos(1));
            end
        end
    else
        %else just use known starting value
        pos_values(1,:) = start_pos;
    end
    
    % get subwindow at current estimated target position, to train classifer
    x = get_subwindow(img, pos, sz, cos_window);
    
    %Kernel Regularized Least-Squares, calculate alphas (in Fourier domain)
    k = dense_gauss_kernel(sigma, x);
    new_alphaf = yf ./ (fft2(k) + lambda);   %(Eq. 7)
    new_z = x;
    
    if i == 1  %first frame, train with a single image
        alphaf = new_alphaf;
        z = x;
    else
        %subsequent frames, interpolate model
        alphaf = (1 - interp_factor) * alphaf + interp_factor * new_alphaf;
        z = (1 - interp_factor) * z + interp_factor * new_z;
    end
    
    %save position
    positions(i,:) = pos;
    
    % visualization
    rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
    if i == 1  %first frame, create GUI
        figure('IntegerHandle','off', 'Name',['Tracker - ' video_path])
        img = rgb2gray(insertText(img,[0 0],i));
        im_handle = imshow(img, 'Border','tight', 'InitialMag',200);
        rect_handle = rectangle('Position',rect_position,'EdgeColor',rect_color,'LineWidth',1);
    else
        try  %subsequent frames, update GUI
            img = rgb2gray(insertText(img,[0 0],i));
            set(im_handle, 'CData', img)
            set(rect_handle, 'Position', rect_position, 'EdgeColor',rect_color,'LineWidth',1);
        catch
            return
        end
    end
    %keyboard;
    drawnow
    pause(0.01)  %uncomment to run slower
    
    if occlude_flag %for screenshot purposes
        pause(2)
        %keyboard;
    end
end

if use_kalman
    kalman_est = positions;
    kalman_psr = PSR_values;
else
    original_est = positions;
    original_psr = PSR_values;
end