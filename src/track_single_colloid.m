%
% Single Particle Tracking Template
%

threshold = 100 ;
frames = 10000 ;

x_cm = zeros(frames,1);
y_cm = zeros(frames,1);

parfor frame = 1 : frames
    
    % Load current image
    img_in = double( imread( [ '/Volumes/usb-gamma/smp_22per_gamma/' num2str( frame ) '.tiff' ] ) ) ;
    
    img_in( img_in < threshold ) = 0 ;
    
    % Coordinate grids
    [X,Y] = meshgrid(1:size(img_in,2), 1:size(img_in,1));

    % Calculate the total mass of the image
    M = sum(img_in(:));

    % Centre of mass
    x_cm(frame) = sum(img_in(:).*X(:)) / M;
    y_cm(frame) = sum(img_in(:).*Y(:)) / M;

end

%%
k_B = 1.380649E-23; T = 300 ; fps = 10;

x_m = (x_cm - mean(x_cm))*0.137e-6;
y_m = (y_cm - mean(y_cm))*0.137e-6;

f = (1:frames)';

trap_arr = [x_m y_m f ones(frames,1)];

[msd_2d, msd_x, msd_y, tau, msd_count] = calcMSD( trap_arr, fps );


[h_x, b_x]=histcounts(x_m);[h_y, b_y]=histcounts(y_m);

h_x_f = imgaussfilt(h_x, 3); h_y_f = imgaussfilt(h_y, 3);

p_x = h_x_f / sum(h_x_f); p_y = h_y_f / sum(h_y_f);

u_x = - log(p_x) * k_B * T; u_y = - log(p_y) * k_B * T;

plot(b_x(1:end-1), u_x,'bo', b_y(1:end-1), u_y, 'ro')

%%
model = @(k,x) k*x.^2;

k0 = 1e-12; % initial guess

k_x = nlinfit(b_x(1:end-1), u_x, model, k0)
k_y = nlinfit(b_y(1:end-1), u_y, model, k0)

%%
x_e = mean(x_m.^2); y_e = mean(y_m.^2);

k_B_x = x_e * k_x / T
