%
% Particle Tracking Tutorial
%
% Origanal: https://site.physics.georgetown.edu/matlab/tutorial.html
%
%

close all ; clear all ;

particles = [] ;
pks = [] ;

for frame = 1 : 100
    
    image_array = double( imread( [ '../tests/data/test_img_' num2str( frame, '%03.f' ) '.tiff' ] ) ) ;

    image_array = image_array( 512 : 612, 512 : 612 ) ;

    excl_dia = 15 ;
    excl_rad = floor( excl_dia / 2 ) ;
    backgrnd = 100 ;


    %   `img_out = bpass( img_in, hpass, lpass, backgrnd, display )`
    filtered_image = bpass_org( image_array, 0, excl_dia, backgrnd ) ;

    % est_pks = pkfnd( img, threshold, excl_dia )
    est_pks = pkfnd_org( image_array, backgrnd, excl_dia ) ;
    pks = [ pks ; est_pks ] ;

    % particles = cntrd( img, est_pks, excl_dia )
    cntrds = cntrd_org( image_array, est_pks, excl_dia ) ;
    cntrds = [ cntrds frame * ones( length( cntrds( : , 1 ) ), 1 ) ] ;
    particles = [ particles ; cntrds ] ;

    % image_array_fig = figure ; colormap('gray'), imagesc( image_array ) ; axis square ;
    % 
    % crc = viscircles( [ cntrds( :, 1 ), cntrds( :, 2 ) ], cntrds( :, 4 ) / 2, 'Color', 'r', 'EnhanceVisibility', false, 'LineWidth', 1 ) ;


end
figure ; plot( particles(:,4), particles(:,3), 'o') 

xyt = particles( :, [ 1 2 5 ] ) ;
param.mem = 4 ;
param.good = 0 ;
param.dim = length( xyt( 1, : ) ) - 1 ;
param.quiet = 0 ;

tic

tracks = track_org( xyt, 13, param ) ;
toc

figure
plot(tracks(:,1),tracks(:,2),'.')

counts = [] ; centers = [] ;
yout = [] ;
xout = [] ;
for n = 1 : max( tracks( :, 4 ) )
    ind = find( tracks( :, 4 ) == n ) ;
    clear counts centers
    dat = tracks(ind, 1 ) - mean(tracks(ind, 1 ) );
    [counts, centers ] = hist( dat, 5 ) ;
    yout = [ yout ; counts ] ;
    xout = [ xout ; centers ] ;
end

plot(rot90(xout(:,1)),rot90(yout(:,1)),'.')


