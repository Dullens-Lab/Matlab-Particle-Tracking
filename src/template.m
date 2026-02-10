
particles = []

parfor i = 1 : 1200

    % Load a single image into the Matlab workspace
    image_raw =  imread(sprintf('/Users/arrancurran/Documents/radboud/Teaching/Soft Matter Practical/Test Data/captured_images_2025-01-25_12-48-01/img_%d.tiff', i) ) ;
    
    % Processed the image using a bpass filter.
    filtered_image = bpass( image_raw, false, true, 140, false) ;
    
    % Display the filtered image
    % imshow( filtered_image )
    
    % Colloid positions, estimated to the nearest pixel
    est_pks = pkfnd( filtered_image, 140, 15) ;
    
    % Find the sub-pixel coordinates for the colloid using a centre-of-mass algorithm.
    cntrds = cntrd( filtered_image, est_pks, 15, true, i ) ;
    
    particles = [ particles ; cntrds ] ;

end
   
    
