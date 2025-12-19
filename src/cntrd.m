%
% Calculates the centroid of bright spots to sub-pixel accuracy.
%
% particles = cntrd( img, est_pks, excl_dia )
% 
% img:          2D array of image pixel values.
%               Particles should be bright spots on a dark background with little noise.
%               Often filtered with bpass() prior to here.
%
% est_pks:      Estimated locations of local maxima to pixel-level accuracy from pkfnd().
%
% excl_dia:     If your data's noisy, (e.g. a single particle has multiple local maxima), then set this optional 
%               keyword to a value slightly larger than the diameter of your blob. If multiple peaks are found 
%               withing a radius of excl_dia/2 then the code will keep only the brightest. Also gets rid of all 
%               peaks within excl_dia of image edge.
%
% returns:      N x 4 array containing, x, y, total brightness and estimated radii for each feature
%
%               particles(:,1) is the x-coordinates.
%               particles(:,2) is the y-coordinates.
%               particles(:,3) is the brightnesses.
%               particles(:,4) is the estimated radii.
%
% NOTES:
%
% If pkfnd(), and cntrd() return more than one location per particle then
% you should try to filter your input more carefully. If you still get
% more than one peak for a particle, use the optional excl_dia parameter 
% in pkfnd().
% 
% If you want sub-pixel accuracy, you need to have a lot of pixels in your 
% window (excl_dia>>1). To check for pixel bias, plot a histogram of the 
% fractional parts of the resulting locations.
%

function cntrds = cntrd( img, est_pks, excl_dia, apply_mask, frame )

    if isa( img, 'double' ) ~= 1, img = double( img ) ; end

    if rem( excl_dia, 2 ) == false 
        warning('Exclusion diameter (excl_dia) must be an odd integer.') ;
        cntrds = [ ] ;
        return ;
    end

    if isempty( est_pks )
        warning('There were no estimated peaks (est_pks) provided. Maybe the threshold in pkfnd() is too high.')
        cntrds = [ ] ;
        return;
    end

    % Number of pixels from the central pixel.
    excl_rad = floor( excl_dia / 2 ) ;

    if apply_mask == true
        % Create a circular mask around estimated peak
        cent_px     = excl_rad + 1 ;
        msk_binary  = zeros( excl_dia ) ;  
        msk_binary( cent_px, cent_px ) = 1 ;    
        msk_binary  = bwdist( msk_binary ) ;
        msk_binary  = msk_binary <= excl_rad ;
    else
        msk_binary  = 1 ;
    end

    msk_ind_x = zeros( excl_dia ) ;

    for n = 1 : excl_dia, msk_ind_x( n, : ) = ( 0 : excl_dia - 1 ) ; end
    
    msk_ind_y = msk_ind_x' ;

    [ est_pks_num, ~ ] = size( est_pks ) ;

    cntrds = zeros( est_pks_num , 5 ) ;

    for n = 1 : est_pks_num

        rows    = est_pks( n, 2 ) - excl_rad  : est_pks( n, 2 ) + excl_rad ;
        cols    = est_pks( n, 1 ) - excl_rad  : est_pks( n, 1 ) + excl_rad ;
        roi     = msk_binary .* img( rows, cols ) ;
        tot_br  = sum( roi, 'all' ) ;

        cntrd_x = est_pks( n, 1 ) + sum( roi .* msk_ind_x, 'all' ) / tot_br - excl_rad ;
        cntrd_y = est_pks( n, 2 ) + sum( roi .* msk_ind_y, 'all' ) / tot_br - excl_rad ;
        pk_val  = max( roi, [], 'all' ) ;
        rad_gyr = 2 * sqrt( sum( roi .^2, 'all' ) / numel( roi ) / pk_val ) ; % Estimate of particles diameter based on radius of gyration
        
        cntrds( n, : ) = [ cntrd_x, cntrd_y, pk_val, rad_gyr, frame ] ;
     
    end