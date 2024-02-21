%
% Particle Tracking Tutorial
%
% Origanal: https://site.physics.georgetown.edu/matlab/tutorial.html
%
%

close all ;

addpath( '../src/', '../img/' )

% image_array = double( imread( 'data/1.tiff' ) ) ;

% image_array = image_array( 1 : 512, 1 : 512 ) ;

% % img = img(1:512,1:512);
% [p,~] = size( image_array ) ;
% for n = 1 : p

%     limit = ( (p / 2 ) - n ) ;
%     pin( :, n ) = linspace( -limit, limit, p ) ;

% end
% pin = abs( pin ) ;
% pin = pin / max(pin, [], 'all' ) ;

% image_array = image_array - (120 * pin) ;

% image_array = image_array - min( image_array, [], 'all' ) ;

% image_array = image_array / max( image_array, [], 'all' ) * 255 ;

excl_dia = 19 ;
excl_rad = floor( excl_dia / 2 ) ;
backgrnd = 120 ;

loop = 1 ;

time_bpass = zeros( loop, 1 ) ;
time_pkfnd = zeros( loop, 1 ) ;
time_cntrd = zeros( loop, 1 ) ;

for n = 1 : loop

    tic
    % `img_out = bpass( img_in, hpass, lpass, backgrnd, display )`
    filtered_image = bpass( image_array, true, 2, backgrnd, false ) ;
    time_bpass( n ) = toc;
    
    image_array_fig = figure ; colormap('gray'), imagesc( image_array ) ; axis square ;
    
    tic
    % est_pks = pkfnd( img, threshold, excl_dia )
    est_pks = pkfnd( filtered_image, backgrnd, excl_dia ) ;
    time_pkfnd = toc ;
    
    tic
    % particles = cntrd( img, est_pks, excl_dia, apply_mask, frame  )
    particles = cntrd( image_array, est_pks, excl_dia, true, 1 ) ;
    time_cntrd( n ) = toc ;
    
end

time_bpass = mean( time_bpass )
time_pkfnd = mean( time_pkfnd )
time_cntrd = mean( time_cntrd )
total = time_bpass + time_pkfnd + time_cntrd
size( particles )

% 
% for p = 1 : length( particles( :, 1 ) )
%     x = particles( p, 1 ) ;
%     y = particles( p, 2 ) ;
%     line( [ x - excl_rad, x + excl_rad ], [ y, y ], 'Color','green' )
%     line( [ x, x ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
%     line( [ x - excl_rad, x - excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
%     line( [ x + excl_rad, x + excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
%     line( [ x - excl_rad, x + excl_rad ], [ y - excl_rad, y - excl_rad ], 'Color','green' )
%     line( [ x + excl_rad, x - excl_rad ], [ y + excl_rad, y + excl_rad ], 'Color','green' )
% end
% 
% crc = viscircles( [ particles( :, 1 ), particles( :, 2 ) ], particles( :, 4 ) / 2, 'Color', 'r', 'EnhanceVisibility', false, 'LineWidth', 1 ) ;

figure ; plot( particles(:,4), particles(:,3), 'o') 