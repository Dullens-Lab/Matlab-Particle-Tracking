%
% Single Particle Tracking Template
%

threshold = 100 ;
frames = 2000 ;

x_cm = zeros(frames,1);
y_cm = zeros(frames,1);

parfor frame = 1 : frames
    
    % Load current image
    img_in = double( imread( [ 'Test Data/captured_images_2025-01-25_12-48-01/img_' num2str( frame ) '.tiff' ] ) ) ;
    
    img_in( img_in < threshold ) = 0 ;
    
    % Coordinate grids
    [X,Y] = meshgrid(1:size(img_in,2), 1:size(img_in,1));

    % Calculate the total mass of the image
    M = sum(img_in(:));

    % Centre of mass
    x_cm(frame) = sum(img_in(:).*X(:)) / M;
    y_cm(frame) = sum(img_in(:).*Y(:)) / M;

end
