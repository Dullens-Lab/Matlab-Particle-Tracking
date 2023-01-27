

clear all ; close all ;
dim = 50 ;
z = 8 ;
[X, Y] = meshgrid( 1 : dim, 1 : dim ) ;
img_gauss = exp( -1 / ( z^2 ) * ( ( Y - dim / 2 ) .^2 + ( X - dim / 2) .^2 ) ) ;


excl_dia = 25 ;


parfor n = 1 : 100000
    noise = abs( randn( size( img_gauss ) ) ) ;
    noise = noise / max( noise, [], 'all' ) ;
    noise = noise / .5 ;
    noise = noise + 1 ;
    

    img_in = img_gauss .* noise ;

    scale2init8 = @( x ) ( x - min( x, [], 'all' ) ) ./ max( ( x - min( x, [], 'all' ) ), [], 'all' ) * 255 ;
    
    img_in  = scale2init8( img_in ) ;

    img_filt = bpass( img_in, false, 8, 50) ;
    img_nfilt = bpass( img_in, false, false, 50) ;

    est_pks_filt = pkfnd( img_filt, false, excl_dia );
    c_filt       = cntrd( img_filt, est_pks_filt, excl_dia, false ) ;

    est_pks_nfilt = pkfnd( img_nfilt, false, excl_dia );
    c_nfilt       = cntrd( img_nfilt, est_pks_nfilt, excl_dia, false ) ;

    x_filt( n ) = c_filt( 1, 1 ) ;
    y_filt( n ) = c_filt( 1, 2 ) ;
    x_nfilt( n ) = c_nfilt( 1, 1 ) ;
    y_nfilt( n ) = c_nfilt( 1, 2) ;

end

[hist_filt, xhist_filt] =  hist( x_filt - mean( x_filt ), 100 ) ;
[hist_nfilt, xhist_nfilt] =  hist( x_nfilt - mean( x_nfilt ), 100 ) ;

% figure ; colormap( 'gray' ) ; imagesc(img_filt) ; hold on; plot(x_filt,y_filt,'o');hold off; 

figure; plot(xhist_filt, hist_filt / max( hist_filt, [], 'all' ),'ro', xhist_nfilt, hist_nfilt / max( hist_nfilt, [], 'all' ),'go')

% d = diff( sort( y ) ) ; 

std(x_filt )
std(x_nfilt ) 

% mean( d( find( d > .1 ) ) ) ;