%
% Finds particle positions in an image to pixel level accuracy. The output here is expected to be passed to cntrd().
%
% After inits, first loops through all peak pixels and checks to see if it is the brightest in a 3 x 3 array (i.e. 8 nearest neighbours).
% Then we exclude all peak pixels whos coordinates lie with the exclusion distance from the image edges.
% Final step is to elimate all but the brightest pixel within the area given by excl_dia.
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
%               keyword to a value slightly larger than the diameter of your blob. If multiple peaks are found 
%               withing a radius of excl_dia/2 then the code will keep only the brightest. Also gets rid of all 
%               peaks within excl_dia of image edge.
%
% returns:      N x 2 array containing, coordinates of local maxima.
%
%               Typically, the return is the input for cntrd().
%
% Inspired by the lmx subroutine of Grier and Crocker's feature.pro
%

%{

NOTES:

When comparing px value with nearest 8 neigbours, the instinct may be to first gather the 8, find the max then compare once, i.e.

    border_pixels = max( [ max( img( row - 1 : 2 : row + 1, range + col ) , [] , 'all' ), max( img( row, col - 1 : 2 : col + 1 ) ) ] ) ;
    if img( row, col ) >= border_pixels

but this is about 5 x slower than just directly comparing every pixel!

CHANGELOG:

Feb 4 2005
Written by Eric R. Dufresne, Yale University.

May 2005
Got rid of ind2rc.m to reduce overhead on tip by Dan Blair. ERD
Added sz keyword. ERD

Jun 2005
Modified to work with one and zero peaks, removed automatic normalization of image. ERD
Due to popular demand, altered output to give x and y instead of row and column. ERD

Aug 24 2005
pkfnd now exits politely if there's nothing above threshold instead of crashing rudely. ERD

Jun 14 2006 
Now exits politely if no maxima found. ERD

Oct 5 2006
Fixed bug that threw away particles with maxima consisting of more than two adjacent points. ERD

Jan 2023
Reformated to meet commenting and nomenclecture standards. AC

%}

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

    if length( pk_px_row ) == 0
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

    % Elimate all but one peak within excl_dia
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

            % if sum( roi( :, 1 ) + roi( :, end ) + roi( 1, : )' + roi( end, : )', 'all' ) ~= 0 
            %     disp('hello from pkfnd')
            % end 

            [ roi_pk_px , roi_pk_px_ind ] = max( roi, [], 'all') ;

            % Here we could ask if the found max is near the centre of the current roi and if not immediatley replace the roi over it.
            % this would prevent ever repetition in the same area, like in the situation where there are enough pk_pxs
            % with a particle such that the first roi only gets some of the pk_pxs and we later will have to repeat.
            
            [ roi_pk_px_row, roi_pk_px_col ] = ind2sub( size( roi ), roi_pk_px_ind ) ;
            
            pk_px_img( rows, cols ) = 0 ;
            
            pk_px_img( pk_px_coords( n, 1 ) - excl_rad + roi_pk_px_row - 1, ...
                       pk_px_coords( n, 2 ) - excl_rad + roi_pk_px_col - 1 )...
                       = roi_pk_px ;

        end

        [ pk_px_row, pk_px_col ] = find( pk_px_img > 0 ) ;

    end

    est_pks = [ pk_px_col, pk_px_row ] ;    