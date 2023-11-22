%
%Three step image manipulation starting with a high frequency pass filter to remove long scale variations. The high pass filtered image, `img_hpass` is then filtered with a low (gaussian) pass filter, to remove pixel noise. Finally, `img_lpass`, has any pixel values below `backgrnd` set to zero. Any step can be skipped with a `false` argument.
%
%   `img_out = bpass( img_in, hpass, lpass, backgrnd, display )`
%
%   `img_in` 2D array of image pixel values.
%
%   `hpass` Set to `true` for highpass filtering. Set to `false` to skip.
%
%   `lpass` Set to `true` to apply a gaussian filter with a strength calculated from the input image. Provide any positive integer for manual control of the gaussian kernel. Set to `false` to skip. For either auto or manual, if the strength of the filter is equal to 1 the image is assumed to be good enough amd gaussian filtering and will be skipped.
%
%   `backgrnd` Reset any pixel values below `backgrnd` to 0. Set to 'false' to skip.
%
%   `display` Plot the image and pixel distribution at each stage of the filtering. Set to `false` or leave blank to skip.
%
%   `img_out` 2D array of filtered image pixel values.
%
%`img_hpass` and `img_lpass` can be returned with `[ img_out, img_hpass ] = bpass()` and `[ img_out, ~, img_lpass ] = bpass()`, respectively.

%{

Notes on convolution:

JWM: Do a 2D convolution with the kernels in two steps each. It is
possible to do the convolution in only one step per kernel with 

  lpass_conv = conv2(lpass_kernel',lpass_kernel,image,'same');
  boxcar_conv = conv2(box_kernel', box_kernel,image,'same');

but for some reason, this is slow. The whole operation could be reduced
to a single step using the associative and distributive properties of
convolution:

  filtered = conv2(image,...
    lpass_kernel'*lpass_kernel - box_kernel'*box_kernel,...
    'same');

But this is also comparatively slow (though inexplicably faster than the
above). It turns out that convolving with a column vector is faster than
convolving with a row vector, so instead of transposing the kernel, the
image is transposed twice.

This is still true as of 2023. Also imgaussfilt() is also slower and has an effect of translating the centroids in X and Y

CHANGELOG:

Feb 1993
Written by David G. Grier, The University of Chicago.

May 1995
Greatly revised version DGG.

Dec 1995
Added /field keyword JCC.

Aug 1999
Memory optimizations and fixed normalization, DGG.

Apr 2004-ish
Converted to Matlab by D.Blair.

June 2005
Fixed some bugs with conv2 to make sure the edges are removed D.B.
Removed inadvertent image shift ERD.

Aug 24 2005
Added threshold to output. Now sets all pixels with negative values equal to zero. Gets rid of ringing 
which was destroying sub-pixel accuracy, unless window size in cntrd was picked perfectly. Now centrd 
gets sub-pixel accuracy much more robustly ERD.

Jun 2007
Refactored for clarity and converted all convolutions to use column vector kernels for speed. Running
on my  macbook, the old version took ~1.3 seconds to do bpass(image,1,19) on a 1024 x 1024 image;
this version takes roughly half that. JWM

Jan 2023
Reformated to meet commenting and nomenclecture standards. AC
Added image scalling on input and ouptut so we process with full bandwidth (0-255) and we return full bandwidth. AC
Added check to see if image in has already been converted to 'double'. AC
Removed edge zero-ing as this happens in pkfnd() by removing peaks within edge rather than blacking out and potentially 
picking parts of particles.

%}

function [ img_out, img_hpass, img_lpass ] = bpass( img_in, hpass, lpass, backgrnd, display )

    if nargin < 4
        warning('No image filtering performed. Not enough arguments provided in bpass( img, hpass, lpass, backgrnd, display )')
        img_out = img_in ;
        return
    end

    if ~exist( 'display', 'var' ), display = false ; end

    if isa( img_in, 'double' ) ~= 1, img_in = double( img_in ) ; end

    normalize   = @( x ) x / sum( x ) ;
    scale2init8 = @( x ) ( x - min( x, [], 'all' ) ) ./ max( ( x - min( x, [], 'all' ) ), [], 'all' ) * 255 ;
    
    img_in      = scale2init8( img_in ) ;
    img_out     = img_in ;
    
        
    %%%     High Pass Filter    %%%
    
    % The kernel is designed to increase the brightness of the center pixel relative to neighboring pixels.
    % The kernel array usually contains a single positive value at its center, which is completely surrounded by negative values.
    % The following array is an example of a 3 by 3 kernel for a high pass filter:
    %
    %   -1/9    -1/9    -1/9
    %   -1/9     8/9    -1/9
    %   -1/9    -1/9    -1/9
    %
    % https://www.l3harrisgeospatial.com/docs/highpassfilter.html
    % TODO: Explore effect of amplitude on kernel
    if hpass
        box_kernel  = - ones( 3 ) / 9 ; box_kernel( 2, 2 ) = 8 / 9 ;
        img_hpass   = conv2( img_out, box_kernel, 'same' ) ;
        img_hpass   = scale2init8( img_hpass ) ;
        img_out     = img_hpass ;
    end

    %%%     Low Pass Filter    %%%
    %
    % The kernel is designed to blur groups of pixels based on lpass. If an integer is not provided, it is estiamted.
    %
    if lpass ~= false

        if islogical( lpass )
            % Estimate noise from input image, see https://doi.org/10.1006/cviu.1996.0060
            % 
            % Noise Estimation Operator, nop
            %    1  -2   1
            %   -2   4  -2
            %    1  -2   1
            %
            nop_bld = [ 1 -2 1 ] ;
            nop     = [ nop_bld ; - nop_bld * 2 ; nop_bld ] ;

            [ img_rows, img_cols ] = size( img_out ) ;
            % TODO: Should we calc nop from the raw input regardless if hpass has happened?
            nop_sigma   = sum( abs( conv2( img_out, nop ) ), 'all'  ) ;
            lpass       = round( nop_sigma * sqrt( .5 * pi ) / ( 6 * ( img_rows - 2 ) * ( img_cols - 2 ) ) ) ;
        end

        if lpass ~= 1 % Dont waste my time with good images!
            lpass_x      = - lpass : lpass ;
            lpass_kernel = normalize( exp( -( lpass_x / ( 2 * lpass ) ) .^2 ) ) ;
            img_lpass    = conv2( img_out, lpass_kernel, 'same' ) ;
            img_lpass    = conv2( img_lpass, lpass_kernel', 'same' ) ;
            img_lpass    = scale2init8( img_lpass ) ;
            img_out      = img_lpass ;
        end
        
    end

    %%%     Zero Background Pixels    %%%

    if backgrnd
        img_base = img_out ;
        img_base( img_base < backgrnd ) = 0 ; 
        img_out = img_base ;
    end


    if display == true

        fov = 36 ;
        figure_img = figure ; colormap( figure_img, 'gray') ; figure_hists = figure ;

        img_hist = @( x )  hist( x, min( x, [], 'all' ) : max( x, [], 'all' ) ) ;

        [ hist_raw, x_hist ] = img_hist( img_in ) ;
        figure_hists ; semilogy1_raw = semilogy( x_hist, sum( hist_raw, 2 ), 'ko' ) ;
        set(semilogy1_raw, 'DisplayName', 'Raw' ) ;
        hold on
        display_raw = subplot( 2, 2, 1, 'Parent', figure_img ) ; image( img_in( 1 : fov, 1 : fov ), 'Parent', display_raw) ;
        title( display_raw, 'Raw Image' ) ; set( display_raw, 'YTickLabel', [ ] ) ; set( display_raw, 'XTickLabel', [ ] ) ;

        if hpass
            [ hist_box, x_hist ] = img_hist( img_hpass ) ;
            figure_hists ; semilogy1_box = semilogy( x_hist, sum( hist_box, 2 ), 'bo' ) ;
            set(semilogy1_box, 'DisplayName', 'Boxcar', 'MarkerFaceColor', 'b' ) ;
            display_box = subplot( 2, 2, 2, 'Parent', figure_img ) ; imagesc( img_hpass( 1 : fov, 1 : fov ), 'Parent', display_box) ;
            title( display_box, 'Boxcar Filtered Image' ) ;set( display_box, 'YTickLabel', [ ] ) ; set( display_box, 'XTickLabel', [ ] ) ;
        end

        if exist( 'img_lpass', 'var' )
            [ hist_g, x_hist ] = img_hist( img_lpass ) ;
            figure_hists ; semilogy1_g = semilogy( x_hist, sum( hist_g, 2 ), 'ro' ) ;
            set( semilogy1_g, 'DisplayName', 'lpassian', 'MarkerFaceColor', 'r' ) ;
            display_lpass = subplot( 2, 2, 3, 'Parent', figure_img ) ; image( img_lpass( 1 : fov, 1 : fov ), 'Parent', display_lpass) ;
            title( display_lpass, 'gaussian Filtered Image' ) ;set( display_lpass, 'YTickLabel', [ ] ) ; set( display_lpass, 'XTickLabel', [ ] ) ;

        end

        [ hist_f, x_hist ] = img_hist( img_out ) ;
        figure_hists ; semilogy1_f = semilogy( x_hist, sum( hist_f, 2 ), 'go' ) ;
        set(semilogy1_f, 'DisplayName', 'Output' ) ;
        display_out = subplot( 2, 2, 4, 'Parent', figure_img ) ; image( img_out( 1 : fov, 1 : fov ), 'Parent', display_out) ;
        title( display_out, 'Output Image' ) ;set( display_out, 'YTickLabel', [ ] ) ; set( display_out, 'XTickLabel', [ ] ) ;

    end