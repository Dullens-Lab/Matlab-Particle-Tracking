%
% Finds particle positions in an image to pixel level accuracy. The output here is expected to be passed to cntrd().
%
% After inits, first loop through all peak pixels and checks to see if it is the brightest in a 3 x 3 array (i.e. 8 nearest neighbours).
% Then exclude all peak pixels who's coordinates lie within the exclusion distance from the image edges.
% Final step is to eliminate all but the brightest pixel within the area given by excl_dia.
%
%
% est_pks = pkfnd( img, threshold, excl_dia )
%
% img:          2D array of image pixel values.
%               Particles should be bright spots on a dark background with little noise.
%               Often filtered with bpass() prior to here.
% 
% threshold:    The minimum brightness of a pixel that might be local maxima.
%               Large values will result in faster code execution but you might miss some particles.
%               Small values will result in slower code execution but you might get some false particles.
%   
% excl_dia:     If your data's noisy, (e.g. a single particle has multiple local maxima), then set this optional 
%               keyword to a value slightly larger than the diameter of your particle. If multiple peaks are found 
%               withing a radius of excl_dia/2 then the code will keep only the brightest. Also gets rid of all 
%               peaks within excl_dia of the image edges.
%
% returns:      N x 2 array containing, coordinates of local maxima.
%
%               Typically, the return is the input for cntrd().
%

function [ est_pks, input_pk_pxs ] = pkfnd( img, th, excl_dia )

    if nargin ~= 3
        warning('Not enough arguemts for pkfnd( img, th, excl_dia )') ;
        est_pks = [ ] ;
        return ;
    end

    if rem( excl_dia, 2 ) == false 
        warning('Exclusion diameter (excl_dia) must be an odd integer.') ;
        est_pks = [ ] ;
        return ;
    end

    if isa( img, 'double' ) ~= 1, img = double( img ) ; end

    [ pk_px_row, pk_px_col ]    = find( img >= th ) ;
    input_pk_pxs                = length( pk_px_row ) ;
    [ img_rows, img_cols ]      = size( img ) ;
    excl_rad                    = floor( excl_dia / 2 ) ;

    if isempty( pk_px_row )
        warning( ['The provided image does not contain any pixel values above the ', num2str(th), ' pixel threshold set in pkfnd()'] ) ;
        est_pks = [ ] ;
        return;
    end

    % Get ride of pks within excl_dia of image edges
    if pk_px_row > 0

        ind = find( pk_px_row > excl_rad & pk_px_row < img_rows - excl_rad ) ;

        pk_px_row = pk_px_row( ind ) ;
        pk_px_col = pk_px_col( ind ) ;

        ind = find( pk_px_col > excl_rad & pk_px_col < img_cols - excl_rad ) ;

        pk_px_col = pk_px_col( ind ) ;
        pk_px_row = pk_px_row( ind ) ;

    end

    pk_px_num       = length( pk_px_row ) ;
    pk_px_coords    = zeros( pk_px_num, 2 ) ;

    cnt = 1 ;

    % Check each pixel above threshold to see if it's brighter than it's 8 neighbors.
    for n = 1 : pk_px_num

        row = pk_px_row( n ) ;
        col = pk_px_col( n ) ;

        if img( row, col ) >= img( row,     col + 1 ) ...
        && img( row, col ) >= img( row - 1, col + 1 ) ...
        && img( row, col ) >= img( row - 1, col     ) ...
        && img( row, col ) >= img( row - 1, col - 1 ) ...
        && img( row, col ) >= img( row,     col - 1 ) ...
        && img( row, col ) >= img( row + 1, col - 1 ) ...
        && img( row, col ) >= img( row + 1, col     ) ...
        && img( row, col ) >= img( row + 1, col + 1 )
        
            pk_px_coords( cnt, : ) = [ row, col ] ; 
            cnt = cnt + 1 ;

        end
    end

    pk_px_coords        = pk_px_coords( 1 : cnt - 1 , : ) ;
    [ pk_px_num, ~ ]    = size( pk_px_coords ) ;

    % Eliminate all but one peak within excl_dia
    if pk_px_num > 1

        pk_px_img = zeros( img_rows, img_cols ) ;
        
        for n = 1 : pk_px_num

            row = pk_px_coords( n, 1 ) ;
            col = pk_px_coords( n, 2 ) ;

            pk_px_img( row, col ) = img( row, col ) ;
        
        end

        for n = 1 : pk_px_num

            rows    = pk_px_coords( n, 1 ) - excl_rad : pk_px_coords( n, 1 ) + excl_rad ;
            cols    = pk_px_coords( n, 2 ) - excl_rad : pk_px_coords( n, 2 ) + excl_rad ;
            roi     = pk_px_img( rows, cols ) ;

            [ roi_pk_px , roi_pk_px_ind ] = max( roi, [], 'all') ;

            [ roi_pk_px_row, roi_pk_px_col ] = ind2sub( size( roi ), roi_pk_px_ind ) ;
            
            pk_px_img( rows, cols ) = 0 ;
            
            pk_px_img( pk_px_coords( n, 1 ) - excl_rad + roi_pk_px_row - 1, ...
                       pk_px_coords( n, 2 ) - excl_rad + roi_pk_px_col - 1 )...
                       = roi_pk_px ;

        end

        [ pk_px_row, pk_px_col ] = find( pk_px_img > 0 ) ;

    end

    est_pks = [ pk_px_col, pk_px_row ] ;    