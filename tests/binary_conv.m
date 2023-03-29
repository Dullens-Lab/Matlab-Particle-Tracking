% Testing the need to conv after binary


addpath( '../src/' ) ;
clear all ; close all ;
scale2init8 = @( x ) ( x - min( x, [], 'all' ) ) ./ max( ( x - min( x, [], 'all' ) ), [], 'all' ) * 255 ;

test_loop   = 100 ;
img_dim     = 40 ;
sigma       = 6 ;
excl_dia    = 31 ;
[X, Y]  = meshgrid( 1 : img_dim, 1 : img_dim ) ;

% Generate fake image
img_in = exp( -1 / ( sigma^2 ) * ( ( Y - img_dim / 2 ) .^2 + ( X - img_dim / 2) .^2 ) ) ;

for p = 1 : test_loop
    noise       = randn( img_dim ) ;
    noise       = noise - min( noise, [], 'all' ) ;
    noise       = noise / max( noise, [], 'all' ) ;
    noise       = noise * 40 + 1 ;

    % Add noise to our image
    img_noise   = img_in .* noise ;
    img_noise   = scale2init8( img_noise ) ;
    
    % Binary image then centroid
    img_binary  = imbinarize( img_noise, 180 );
    est_pks_b   = pkfnd( img_binary, 0, excl_dia );
    cntrd_b     = cntrd( img_binary, est_pks_b, excl_dia, false ) ;
    x_b( p )    = cntrd_b( 1, 1 ) ;
    y_b( p )    = cntrd_b( 1, 2 ) ;

    % Binary then lpass image then centroid
    img_lpass   = bpass( img_binary, false, 4, false) ;
    est_pks_f   = pkfnd( img_lpass, 0, excl_dia );
    cntrd_f     = cntrd( img_lpass, est_pks_f, excl_dia, false ) ;
    x_f( p )    = cntrd_f( 1, 1 ) ;
    y_f( p )    = cntrd_f( 1, 2 ) ;

end


mean_x = mean( x_b ./ x_f ) ;
mean_y = mean( y_b ./ y_f ) ;

if mean_x == 1 && mean_y == 1
    disp('The two methods are identical' )
else
    disp('Someting is different')
end

