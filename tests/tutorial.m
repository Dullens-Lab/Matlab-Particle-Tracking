%
% Particle Tracking Tutorial
%
% Origanal: https://site.physics.georgetown.edu/matlab/tutorial.html
%
%

close all ; clear all ;

addpath( '../src/' )

particles = [] ;

for frame = 1 : 2
    
    image_array = double( imread( [ 'data/' num2str( frame ) '.tiff' ] ) ) ;

    excl_dia = 15 ;
    excl_rad = floor( excl_dia / 2 ) ;
    backgrnd = 120 ;

    %   `img_out = bpass( img_in, hpass, lpass, backgrnd, display )`
    filtered_image = bpass( image_array, false, true, backgrnd, true ) ;

    % est_pks = pkfnd( img, threshold, excl_dia )
    est_pks = pkfnd( filtered_image, backgrnd, excl_dia ) ;

    % particles = cntrd( img, est_pks, excl_dia )

    cntrds = cntrd( image_array, est_pks, excl_dia, true, frame ) ;
    particles = [ particles ; cntrds ] ;

    image_array_fig = figure ; colormap('gray'), imagesc( image_array ) ; axis square ;

    crc = viscircles( [ cntrds( :, 1 ), cntrds( :, 2 ) ], cntrds( :, 4 ) / 2, 'Color', 'r', 'EnhanceVisibility', false, 'LineWidth', 1 ) ;


end

figure ; plot( particles(:,4), particles(:,3), 'o') 

xyt = particles( :, [ 1 2 5 ] ) ;

tracks = track( xyt, 5 ) ;
