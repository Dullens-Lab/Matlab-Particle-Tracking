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

    image_array = image_array( 512 : 612, 512 : 612 ) ;

    % img_out = bpass( img_in, hpass, lpass, backgrnd, display )
    filtered_image = bpass_org( image_array, 0, excl_dia, backgrnd ) ;

    % est_pks = pkfnd( img, threshold, excl_dia )
    est_pks = pkfnd_org( filtered_image, backgrnd, excl_dia ) ;
    pks = [ pks ; est_pks ] ;

    % particles = cntrd( img, est_pks, excl_dia )
    cntrds = cntrd_org( filtered_image, est_pks, excl_dia ) ;
    cntrds = [ cntrds frame * ones( length( cntrds( : , 1 ) ), 1 ) ] ;
    particles = [ particles ; cntrds ] ;

end

xyt = particles( :, [ 1 2 5 ] ) ;

tracks = track_org( xyt, maxdisp, param ) ;

counts = [] ; centers = [] ;
yout = [] ; xout = [] ;

for n = 1 : max( tracks( :, 4 ) )
    ind = find( tracks( :, 4 ) == n ) ;
    clear counts centers
    dat = tracks(ind, 1 ) - mean(tracks(ind, 1 ) );
    [counts, centers ] = hist( dat, 5 ) ;
    yout = [ yout ; counts ] ;
    xout = [ xout ; centers ] ;
end

figure ; plot( particles(:,4), particles(:,3), 'o') 
figure ; plot(tracks(:,1),tracks(:,2),'.')
figure ; plot(rot90(xout(:,1)),rot90(yout(:,1)),'.')


