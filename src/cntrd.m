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

%{

CHANGELOG:

Feb 4 2005
Written by Eric R. Dufresne, Yale University.

May 2005
Inputs diamter instead of radius.

Jun 2005
Added code from imdist/dist to make this stand alone. DB

Increased frame of reject locations around edge to 1.5*excl_dia ERD

By popular demand, 
1. altered input to be formatted in x,y space instead of row, column space. ERD
2. added forth column of output, rg^2. ERD

Aug 2005
Outputs had been shifted by [0.5,0.5] pixels.  No more! ERD

Aug 24 2005
Woops!  That last one was a red herring.  The real problem is the "ringing" from the output of bpass.
I fixed bpass (see note), and no longer need this kludge. Also, made it quite nice if est_pks=[]; ERD

Jun 2006
Added size and brightness output ot interactive mode. Also fixed bug in calculation of rg^2. ERD

Jun 2007
Small corrections to documentation. JWM

Jan 2023
Reformated to meet commenting and nomenclecture standards. AC
Removed interactive option. AC
Changed excl_rad such that we consider n pixels from the given peak where 2n + 1 is the input excl_dia.
Removed filtering image edges of peaks and this is already happening in pkfnd(). AC
Changed radius of giration calc. I dont trust rg = ( sum( roi .* dst2, 'all' ) / norm ) ; since the applied mask influences
the estimated radius.
Changed the mask such that it is a pixellated circle of diameter = excl_dia rather than excl_dia - 1.
Added option to include mask or not.



    % Create mask - window around trial location over which to calculate the centroid

    % msk_range = ( - excl_rad : excl_rad ) .^2 ;
    
    % cent_px = excl_rad + 1 ;

    % msk_inv = zeros( excl_dia ) ;
    
    % for i = 1 : excl_dia
    %     msk_inv( i, : ) = sqrt( ( i - cent_px ) ^2 + msk_range ) ;
    % end

    % ind = find( msk_inv <= excl_rad ) ;

    % msk_binary = zeros( excl_dia ) ;

    % msk_binary( ind ) = 1.0 ;

    % dst2 = msk_binary .* ( msk_inv .^2 ) ;

%}
function cntrds = cntrd( img, est_pks, excl_dia, apply_mask )

    if isa( img, 'double' ) ~= 1, img = double( img ) ; end

    if rem( excl_dia, 2 ) == false 
        warning('Exclusion diameter (excl_dia) must be an odd integer.') ;
        particles = [ ] ;
        return ;
    end

    if isempty( est_pks )
        warning('There were no estimated peaks (est_pks) provided. Maybe the threshold in pkfnd() is too high.')
        particles = [ ] ;
        return;
    end

    % Number of pixels from the central pixel.
    excl_rad = floor( excl_dia / 2 ) ;

    if apply_mask == true
        % Create a circular mask around trial
        cent_px     = excl_rad + 1 ;
        msk_binary  = zeros( excl_dia ) ;  
        msk_binary( cent_px, cent_px ) = 1 ;    
        msk_binary  = bwdist( msk_binary );
        msk_binary  = msk_binary <= excl_rad;
    else
        msk_binary  = 1 ;
    end

    msk_ind_x = zeros( excl_dia ) ;

    for n = 1 : excl_dia, msk_ind_x( n, : ) = ( 0 : excl_dia - 1 ) ; end
    
    msk_ind_y = msk_ind_x' ;

    [ est_pks_num, ~ ] = size( est_pks ) ;

    cntrds = zeros( est_pks_num , 4 ) ;

    for n = 1 : est_pks_num
        rows    = est_pks( n, 2 ) - excl_rad  : est_pks( n, 2 ) + excl_rad ;
        cols    = est_pks( n, 1 ) - excl_rad  : est_pks( n, 1 ) + excl_rad ;
        roi     = msk_binary .* img( rows, cols ) ;
        tot_br  = sum( roi, 'all' ) ;

        cntrd_x = est_pks( n, 1 ) + sum( roi .* msk_ind_x, 'all' ) / tot_br - excl_rad ;
        cntrd_y = est_pks( n, 2 ) + sum( roi .* msk_ind_y, 'all' ) / tot_br - excl_rad ;
        pk_val  = max( roi, [], 'all' ) ;
        rad_gyr = 2 * sqrt( sum( roi .^2, 'all' ) / numel( roi ) / pk_val ) ; % Estimate of particles diameter based on radius of gyration
        
        cntrds( n, : ) = [ cntrd_x, cntrd_y, pk_val, rad_gyr ] ;
     
    end