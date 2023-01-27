

%img_trim = img_in( 50:85, 50:85);

dim = 100 ;
z = 20 ;
[X, Y] = meshgrid(1:dim, 1:dim) ;
img_trim = exp(-1/(z^2)*((Y-dim/2).^2 + (X-dim/2).^2));


for n = 1 : 10
    pnoise = rand(1) - .5 ;
    [X, Y] = meshgrid(1:dim, (1:dim)+10*pnoise) ;
    img_trim = exp(-1/(z^2)*((Y-dim/2).^2 + (X-dim/2).^2));

    %noise = .05 * rand( size( img_trim ) );

    img_in_noise = img_trim;% + noise ;

    scale2init8 = @( x ) ( x - min( x, [], 'all' ) ) ./ max( ( x - min( x, [], 'all' ) ), [], 'all' ) * 255 ;
    img_in_noise         = scale2init8( img_in_noise ) ;

    img_out = img_in_noise ; %bpass( img_in_noise, 0) ; 
    est_pks = pkfnd(img_out,0, 41 );

    c = cntrd( img_in_noise, est_pks, 41, false ) ;
    %cntrds{ n } = c ;
    pn(n) = pnoise ;
    x(n) = c(1,1);
    y(n) = c(1,2);

end

hist( y, 100 ) ;

figure; imagesc(img_in_noise) ; hold on; plot(x,y,'o');hold off

%figure;plot(x,y,'o')