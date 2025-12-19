
COORD = struct ;
i = 1 ;

% Load a single image into the Matlab workspace
image_raw =  imread(' ') ;

% Processed the image using a bpass filter.
filtered_image = bpass( image_raw, true, true, 0, false) ;

% Display the filtered image
imshow( filtered_image )

% Colloid positions, estimated to the nearest pixel
est_pks = pkfnd( filtered_image, 10, 3 ) ;

% Find the sub-pixel coordinates for the colloid using a centre-of-mass algorithm.
cntrds = cntrd( filtered_image, est_pks, 5, true, 1 ) ;

% Append frame number to the coordinates.
cntrds = [ cntrds i * ones( size( cntrds, 1 ), 1 ) ] ;

% Store the centroids in a structure.
COORD( i ).cntrds = cntrds ;
