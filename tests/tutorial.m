%
% Particle Tracking Tutorial
%
% Origanal: https://site.physics.georgetown.edu/matlab/tutorial.html
%
%

close all ; clear all ;

addpath( '../src/' )

particles = [] ;
pks = [] ;

for frame = 1 : 100
    
    image_array = double( imread( [ 'data/' num2str( frame ) '.tiff' ] ) ) ;

    excl_dia = 15 ;
    excl_rad = floor( excl_dia / 2 ) ;
    backgrnd = 120 ;


    %   `img_out = bpass( img_in, hpass, lpass, backgrnd, display )`
    filtered_image = bpass( image_array, false, false, backgrnd, false ) ;

    % est_pks = pkfnd( img, threshold, excl_dia )
    est_pks = pkfnd( image_array, backgrnd, excl_dia ) ;
    pks = [ pks ; est_pks ] ;

    % particles = cntrd( img, est_pks, excl_dia )
    cntrds = cntrd( image_array, est_pks, excl_dia, false, frame ) ;
    particles = [ particles ; cntrds ] ;

    % image_array_fig = figure ; colormap('gray'), imagesc( image_array ) ; axis square ;
    % 
    % crc = viscircles( [ cntrds( :, 1 ), cntrds( :, 2 ) ], cntrds( :, 4 ) / 2, 'Color', 'r', 'EnhanceVisibility', false, 'LineWidth', 1 ) ;


end

figure ; plot( particles(:,4), particles(:,3), 'o') 

xyt = particles( :, [ 1 2 5 ] ) ;
tic
tracks = track( xyt, 5 ) ;
toc

figure
plot(tracks(:,1),tracks(:,2),'.')