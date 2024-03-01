%
% Particle Tracking Tutorial
%
% Origanal: https://site.physics.georgetown.edu/matlab/tutorial.html
%
%

close all ; clearvars ;

addpath( '../' )
tutorial_params

for frame = 1 : 100
    
    image_array = double( imread( [ '../tests/data/test_img_' num2str( frame, '%03.f' ) '.tiff' ] ) ) ;

    image_array = image_array( 412 : 612, 412 : 612 ) ;

    % img_out = bpass( img_in, hpass, lpass, backgrnd, display )
    filtered_image = bpass( image_array, false, false, backgrnd, false ) ;

    % est_pks = pkfnd( img, threshold, excl_dia )
    est_pks = pkfnd( image_array, backgrnd, excl_dia ) ;
    pks = [ pks ; est_pks ] ;

    % particles = cntrd( img, est_pks, excl_dia )
    cntrds = cntrd( image_array, est_pks, excl_dia, true, frame ) ;
    particles = [ particles ; cntrds ] ;

end

xyt = particles( :, [ 1 2 5 ] ) ;

tracks = track( xyt, 11, param ) ;

counts = [] ; centers = [] ;
yout = [] ; xout = [] ;

for n = 1 : max( tracks( :, 4 ) )
    ind = find( tracks( :, 4 ) == n ) ;
    clear counts centers
    dat = tracks(ind, 1 ) - mean(tracks(ind, 1 ) );
    [counts, centers ] = hist( dat, 20 ) ;
    yout = [ yout ; counts ] ;
    xout = [ xout ; centers ] ;
end

figure ; plot( particles(:,4), particles(:,3), 'o') 
figure ; plot(tracks(:,1),tracks(:,2),'.')

ind = 1 : max( tracks( :, 4 ) ) ;
figure ; plot(mean( xout( ind, : ) ), mean( yout( ind, : ) ), 'o')



