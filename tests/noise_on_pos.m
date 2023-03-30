
addpath( '../src/' ) ;
clear all ; close all ;
dim     = 30 ;
sigma   = 4 ;
[X, Y]  = meshgrid( 1 : dim, 1 : dim ) ;

img_in = exp( -1 / ( sigma^2 ) * ( ( Y - dim / 2 ) .^2 + ( X - dim / 2) .^2 ) ) ;


excl_dia = 29 ;
excl_rad = floor( excl_dia / 2 ) ;

scale2init8 = @( x ) ( x - min( x, [], 'all' ) ) ./ max( ( x - min( x, [], 'all' ) ), [], 'all' ) * 255 ;

img_in  = scale2init8( img_in ) ;

img_in  = double( uint8( img_in ) ) ;

for t = 1 : 1

for n = 1 : 10000

    noise = rand( dim ) ;
    noise = noise - min( noise, [], 'all' ) ;
    noise = noise / max( noise, [], 'all' ) ;
    noise = noise * 40 ;
    
        
    img_noise = img_in + noise ;
    img_noise   = scale2init8( img_noise ) ;

    img_filt    = bpass( img_noise, false, 5, 80) ;
    img_filt    = uint8( img_filt );

    img_nfilt   = bpass( img_noise, false, false, 40) ;
    img_nfilt    = uint8( img_nfilt );

    est_pks_f   = pkfnd( img_filt, 0, excl_dia );
    cntrd_f     = cntrd( img_filt, est_pks_f, excl_dia, false ) ;

    est_pks_nf  = pkfnd( img_nfilt, 0, excl_dia );
    cntrd_nf    = cntrd( img_nfilt, est_pks_nf, excl_dia, false ) ;

    x_f( t, n )    = cntrd_f( 1, 1 ) ;
    y_f( t, n )    = cntrd_f( 1, 2 ) ;

    x_nf( t, n )   = cntrd_nf( 1, 1 ) ;
    y_nf( t, n )   = cntrd_nf( 1, 2 ) ;

end

x_f_mean = x_f - mean( x_f ) ;
y_f_mean = y_f - mean( y_f ) ;

x_nf_mean = x_nf - mean( x_nf ) ;
y_nf_mean = y_nf - mean( y_nf ) ;


var_f( t ) = std( x_f ) ;

var_nf( t ) = std( x_nf ) ; 

end

[hist_f, xhist_f]   =  hist( x_f, 100 ) ;
[hist_nf, xhist_nf] =  hist( x_nf, 100 ) ;

figure ; colormap( 'gray' ) ; imagesc( img_nfilt ) ; hold on ; plot( x_nf, y_nf, 'o' ) ;
x = dim / 2 ;
y = x ;
line( [ x - excl_rad, x + excl_rad ], [ y, y ], 'Color','green' )
line( [ x, x ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
line( [ x - excl_rad, x - excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
line( [ x + excl_rad, x + excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
line( [ x - excl_rad, x + excl_rad ], [ y - excl_rad, y - excl_rad ], 'Color','green' )
line( [ x + excl_rad, x - excl_rad ], [ y + excl_rad, y + excl_rad ], 'Color','green' )
hold off ; 

figure ; colormap( 'gray' ) ; imagesc( img_filt ) ; hold on ; plot( x_f, y_f, 'o' ) ;
% x = mean( x_f ) ;
% y = mean( y_f ) ;
line( [ x - excl_rad, x + excl_rad ], [ y, y ], 'Color','green' )
line( [ x, x ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
line( [ x - excl_rad, x - excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
line( [ x + excl_rad, x + excl_rad ], [ y - excl_rad, y + excl_rad ], 'Color','green' )
line( [ x - excl_rad, x + excl_rad ], [ y - excl_rad, y - excl_rad ], 'Color','green' )
line( [ x + excl_rad, x - excl_rad ], [ y + excl_rad, y + excl_rad ], 'Color','green' )
hold off ; 

figure ; plot( xhist_f, hist_f / max( hist_f, [], 'all' ), 'ro', xhist_nf, hist_nf / max( hist_nf, [], 'all' ), 'go' )

disp( [ num2str( std( x_nf ) ), ' - Unfiltered Variance' ] )
disp( [ num2str( mean( x_nf ) ), ' - Unfiltered Mean' ] )

disp( [ num2str( std( x_f ) ), ' - Filtered Variance ' ] )
disp( [ num2str( mean( x_f ) ), ' - Filtered Mean ' ] )

img_out_noise = imresize( img_noise,  [ 100 100 ] ) ;
img_out_filt = imresize( img_filt,  [ 100 100 ] ) ;

imwrite( uint8( img_out_noise ), 'img/img_in_noisy.jpg', 'jpeg' ) ;
imwrite( uint8( img_out_filt ), 'img/img_in_noisy_filtered.jpg', 'jpeg' ) ;