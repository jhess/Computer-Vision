
%  Exploiting the Circulant Structure of Tracking-by-detection with Kernels
%
%  Main script for tracking, with a gaussian kernel.
%
%  Jo?o F. Henriques, 2012
%  http://www.isr.uc.pt/~henriques/

clc; clear; warning off; close all
%choose the path to the videos (you'll be able to choose one with the GUI)
base_path = './data/';

% occlusion detection parameters
peak_sz = [21, 21];
occlusion_threshold = 4.6;
diagnostic_plots = false;

%parameters according to the paper
padding = 1;					%extra area surrounding the target
output_sigma_factor = 1/16;		%spatial bandwidth (proportional to target)
sigma = 0.2;					%gaussian kernel bandwidth
lambda = 1e-2;					%regularization
interp_factor = 0.075;			%linear interpolation factor for adaptation

%pos = box_position(1:2)+floor(box_position(3:4)/2);
%harcoded position of box for tracking
start_pos = [65 135];




%notation: variables ending with f are in the frequency domain.

%ask the user for the video
video_path = choose_video(base_path);
if isempty(video_path), return, end  %user cancelled
[img_files, pos, target_sz, resize_image, ground_truth, video_path] = ...
	load_video_info(video_path);

target_sz = [42 38];

%window size, taking padding into account
sz = floor(target_sz * (1 + padding));

%desired output (gaussian shaped), bandwidth proportional to target size
output_sigma = sqrt(prod(target_sz)) * output_sigma_factor;
[rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
y = exp(-0.5 / output_sigma^2 * (rs.^2 + cs.^2));
yf = fft2(y);

%store pre-computed cosine window
cos_window = hann(sz(1)) * hann(sz(2))';

pos = start_pos;

time = 0;  %to calculate FPS
positions = zeros(numel(img_files), 2);  %to calculate precision

num_frames = numel(img_files);
psr_values = zeros(num_frames, 1);
pos_values = [];
% array of occlusion frames when using kalman filter
occlusion_frames = [];

for frame = 1:numel(img_files)
    %% load and process image
	% load image
	im = imread([video_path img_files{frame}]);
    % convert to grayscale
    if size(im,3) > 1
		im = rgb2gray(im);
    end
    % resize
	if resize_image
		im = imresize(im, 0.5);
	end
	
	tic()
	
	%% extract and pre-process subwindow
	x = get_subwindow(im, pos, sz, cos_window);
	
    rect_color = 'g';
    
    %% locate target in the subwindow
	if frame > 1
		% calculate response of the classifier at all locations
		k = dense_gauss_kernel(sigma, x, z);
		response = real(ifft2(alphaf .* fft2(k)));   %(Eq. 9)
        
        % target location is at the maximum response
		[row, col] = find(response == max(response(:)), 1);
		pos = pos - floor(sz/2) + [row, col];
        pos_values(frame,:) = pos;
        
        % check for occlusion 
        %[gmax, gmax_idx] = max(response(:));
        gmax = response(row,col);
        
        %[row, col] = ind2sub(size(response), gmax_idx);
%         peak_c1 = [row - (peak_sz(1) - 1)/2, col - (peak_sz(2) - 1)/2];
%         peak_c2 = peak_c1 + peak_sz;
%         sidelobe_region = true(size(response));
%         sidelobe_region(peak_c1(1):peak_c2(1), peak_c1(2):peak_c2(2)) = false;
%         sidelobe_mean = mean(response(sidelobe_region));
%         sidelobe_stdev = std(response(sidelobe_region));
            window_size = round(min(sz)/7);

            %sidelobe region target area - peak window area
            sidelobe_region = zeros(1,numel(response)-window_size^2);
            index = 1;
            sq_radius = (window_size-1)/2;
            %iterate for each row and column of the response to calculate
            %the sidelobe regions
            for i = 1:size(response,1)
                for j = 1:size(response,2)
                    %sidelobe consists of values outside target window around max response
                    if (i < row - sq_radius || i > row + sq_radius || j < col - sq_radius || j > col + sq_radius)
                        sidelobe_region(index) = response(i,j);
                        index = index+1;
                    end
                end
            end
            sidelobe_mean = mean(sidelobe_region);
            sidelobe_stdev = std(sidelobe_region);
        
        % calculate psr
        psr = (gmax - sidelobe_mean) / sidelobe_stdev;
        psr_values(frame) = psr;
        
        if diagnostic_plots
            figure;
            imdisp(response);
            hold on;
            rectangle('Position', [peak_c1([2, 1]), peak_sz([2, 1])], 'EdgeColor', 'g');
        end
        
        %check for occlusion based on psr
        if psr < occlusion_threshold
                occlusion_frames = [occlusion_frames frame];
                occlude_flag = true;
                rect_color = 'r';
                %pause(0.5)
                fprintf('Occlusion detected in frame %d, psr = %g\n',frame,psr)
                
                %kalman filter prediction
                %initialize num of samples
                samples = 10;
                
                %frame
                f = 1; 
                t = 0 : f : f*samples;
                
                %observed values for x and y
                observed_x = (pos_values(end-1,2) - pos_values(end-samples,2))/(samples-1);
                observed_y = (pos_values(end-1,1) - pos_values(end-samples,1))/(samples-1);
                
                %initial values
                x_init = pos_values(end-samples,2);
                y_init = pos_values(end-samples,1);
                
                %velocity values for x and y
                x_true = x_init + observed_x * t;
                y_true = y_init + observed_y * t;
                
                %kalman x and y values and previous values
                Xk_prev = [0;.5*observed_x];
                Yk_prev = [0;.5*observed_y];
                Xk = [];
                Yk = [];
                
                scale = [1 f;0  1];
                s = 1;
                % G1 is the recursive gain matrix used later 
                % during the kalman iteration
                G1 = [s^2 0;0 s^2];
                
                % Q is the process noise covariance matrix
                Q = [0 0; 0 0];
                
                % H is measurement matrix
                H = [1 0];
                
                sigma_meas = 1;
                R = sigma_meas^2;
                
                % buffers for kalman filter iteration
                Xk_2 = zeros(2,samples+1);
                Yk_2 = zeros(2,samples+1);
                
                Xk_2(:,1) = Xk_prev;
                Yk_2(:,1) = Yk_prev;
                
                Y_2_x = zeros(1,samples+1);
                Y_2_y = zeros(1,samples+1);
                
                for h=1:samples
                    
                    % Y is the measurement vector
                    % Y is equal to observed data + uncertainty from occlusion      
                    Y_x = x_true(h+1)+sigma_meas*randn;
                    Y_y = y_true(h+1)+sigma_meas*randn;
                    %updated measurement vector values
                    Y_2_x(h+1) = Y_x;
                    Y_2_y(h+1) = Y_y;
                    
                    % Kalman iteration covariance matrix
                    P = scale*G1*scale' + Q;
                    
                    %Sub equation for inverse = H*P1*H' + R;
                    
                    % G is the gain matrix 
                    G = P*H'*inv(H*P*H' + R);
                    %recursive Gain matrix
                    G1 = P - G*H*P;
                    
                    Xk = scale*Xk_prev + G*(Y_x-H*scale*Xk_prev);
                    Yk = scale*Yk_prev + G*(Y_y-H*scale*Yk_prev);
                    %updated measurement vector values using sample
                    %iteration
                    Xk_2(:,h+1) = Xk;
                    Yk_2(:,h+1) = Yk;
                    
                    % store measurement vector values for the next iteration
                    Xk_prev = Xk;
                    Yk_prev = Yk;
                end
                %store the position of the measurement vector for the frame
                pos = [Yk_2(1,end) Xk_2(1,end)];
                %update position with kalman filter estimate
                pos_values(frame,2) = round(pos(2));
                pos_values(frame,1) = round(pos(1));
        end
    else
        %if first frame just use known starting value
        pos_values(1,:) = start_pos;
    end
	
	%% get subwindow at current estimated target position, to train classifer
	x = get_subwindow(im, pos, sz, cos_window);
	
	%% Kernel Regularized Least-Squares, calculate alphas (in Fourier domain)
	k = dense_gauss_kernel(sigma, x);
	new_alphaf = yf ./ (fft2(k) + lambda);   %(Eq. 7)
	new_z = x;
	
	if frame == 1  %first frame, train with a single image
		alphaf = new_alphaf;
		z = x;
	else
		%subsequent frames, interpolate model
		alphaf = (1 - interp_factor) * alphaf + interp_factor * new_alphaf;
		z = (1 - interp_factor) * z + interp_factor * new_z;
	end
	
	%% save position and calculate FPS
	positions(frame,:) = pos;
	time = time + toc();
	
	%% visualization
	rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
    if frame == 1  %first frame, create GUI
        figure('IntegerHandle','off', 'Name',['Tracker - ' video_path])
        im = rgb2gray(insertText(im,[0 0],frame));
        im_handle = imshow(im, 'Border','tight', 'InitialMag',200);
        rect_handle = rectangle('Position',rect_position,'EdgeColor',rect_color,'LineWidth',1);
    else
        try  %subsequent frames, update GUI
            im = rgb2gray(insertText(im,[0 0],frame));
            set(im_handle, 'CData', im)
            set(rect_handle, 'Position', rect_position, 'EdgeColor',rect_color,'LineWidth',1);
        catch
            return
        end
    end
% 	if frame == 1  %first frame, create GUI
% 		figure('Name',['Tracker - ' video_path])
%         subplot(2, 1, 1);
% 		im_handle = imshow(im, 'Border','tight', 'InitialMag',200);
% 		rect_handle = rectangle('Position',rect_position, 'EdgeColor','g');
%         subplot(2, 1, 2);
%         psr_plot_handle = plot(psr_values);
%         xlabel('Frame');
%         ylabel('psr');
% 	else
% 		try  %subsequent frames, update GUI
% 			set(im_handle, 'CData', im)
% 			set(rect_handle, 'Position', rect_position)
%             set(psr_plot_handle, 'YData', psr_values);
% 		catch  %, user has closed the window
% 			return
% 		end
% 	end
	
	drawnow
% 	pause(0.05)  %uncomment to run slower
end

%store kalman filter psr positions and values
kalman_est = positions;
kalman_psr = psr_values;

if resize_image, positions = positions * 2; end

disp(['Frames-per-second: ' num2str(numel(img_files) / time)])

%show the precisions plot
show_precision(positions, ground_truth, video_path)

figure;
plot(psr_values);
xlabel('Frame');
ylabel('psr');
