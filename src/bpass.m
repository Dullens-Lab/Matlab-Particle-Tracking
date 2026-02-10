%
% Three step image manipulation starting with a high frequency pass filter to remove long scale variations. The high pass filtered image, `img_hpass` is then filtered with a low (gaussian) pass filter, to remove pixel noise. Finally, `img_lpass`, has any pixel values below `backgrnd` set to zero. Any step can be skipped with a `false` argument.
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
%   `img_hpass` and `img_lpass` can be returned with `[ img_out, img_hpass ] = bpass()` and `[ img_out, ~, img_lpass ] = bpass()`, respectively.

function img_out = bpass( img_in, lpass, backgrnd, display )
    
    sub_plots = 2 ;

    if nargin < 3
        warning('No image filtering performed. Not enough arguments provided in bpass( img, hpass, lpass, backgrnd, display )')
        img_out = img_in ;
        return
    end

    if ~exist( 'display', 'var' ), display = false ; end
    % Convert to double 
    if isa( img_in, 'double' ) ~= 1, img_in = double( img_in ) ; end

    normalize   = @( x ) x / sum( x ) ;
    scale2init8 = @( x ) x ; %( x - min( x, [], 'all' ) ) ./ max( ( x - min( x, [], 'all' ) ), [], 'all' ) * 255 ;
    
    % NOTE: This can be problematic in the scenario where we have a group
    % of images but some images contain no objects. In this case we scale
    % the noise to 255 and therefore pick up these bright noisy pixels

    img_in      = scale2init8( img_in ) ;
    img_out     = img_in ;
    
    %%%     Low Pass Filter    %%%
    %
    % The kernel is designed to blur groups of pixels based on lpass. If an integer is not provided, it is estimated.
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
        
        sub_plots = 3 ;

    end

    %%%     Zero Background Pixels    %%%

    if backgrnd
        img_base = img_out ;
        img_base( img_base < backgrnd ) = 0 ;
        img_base = scale2init8( img_base ) ;
        img_out = img_base ;
    end

    if display == true

        fov = 512 ;
        figure_img = figure ; colormap( figure_img, 'gray') ; axis equal ; figure_hists = figure ;
        make_square = @(ax) set(ax, 'DataAspectRatio',[1 1 1], ...   % square pixels
                           'PlotBoxAspectRatio',[1 1 1]);    % square axes box

        img_hist = @( x )  hist( x, min( x, [], 'all' ) : max( x, [], 'all' ) ) ;
        
        axh = axes('Parent', figure_hists);  % histogram axes
        hold(axh,'on')
        set(axh, 'YScale', 'log');
        set(axh, 'Box','on', ...
                 'TickDir','in', ...
                 'XMinorTick','on', ...
                 'YMinorTick','on');
        
        axh.XRuler.TickLabelGapOffset = 0;  % harmless; keeps layout sane
        axh.YRuler.TickLabelGapOffset = 0;
        
        h = gobjects(0);   % line handles


        [ hist_raw, x_hist ] = img_hist( img_in ) ;
        h(end+1) = semilogy(axh, x_hist, sum(hist_raw,2), 'ko',...
            'DisplayName','Input Image', 'MarkerFaceColor','k');

        display_raw = subplot( 1, sub_plots, 1, 'Parent', figure_img ) ; image( img_in( 1 : fov, 1 : fov ), 'Parent', display_raw) ; make_square(display_raw) ;
        title( display_raw, 'Raw Image' ) ; set( display_raw, 'YTickLabel', [ ] ) ; set( display_raw, 'XTickLabel', [ ] ) ;

        if exist( 'img_lpass', 'var' )
            [ hist_g, x_hist ] = img_hist( img_lpass ) ;
            h(end+1) = semilogy(axh, x_hist, sum( hist_g, 2 ), 'ro',...
                        'DisplayName','Gaussian (Low Pass)', 'MarkerFaceColor','r');

            display_lpass = subplot( 1, sub_plots, 2, 'Parent', figure_img ) ; image( img_lpass( 1 : fov, 1 : fov ), 'Parent', display_lpass) ; make_square(display_lpass) ;
            title( display_lpass, 'Gaussian (Low Pass) Filtered Image' ) ;set( display_lpass, 'YTickLabel', [ ] ) ; set( display_lpass, 'XTickLabel', [ ] ) ;
        end

        [ hist_f, x_hist ] = img_hist( img_out ) ;
        h(end+1) = semilogy(axh, x_hist, sum(hist_f,2), 'go',...
            'DisplayName','Output', 'MarkerFaceColor','g');

        xlabel('Pixel Value', 'Interpreter', 'latex')
        ylabel('Pixel Count', 'Interpreter', 'latex')
        legend(axh, h, 'Location','best');

        display_out = subplot( 1, sub_plots, sub_plots, 'Parent', figure_img ) ; image( img_out( 1 : fov, 1 : fov ), 'Parent', display_out) ; make_square(display_out) ;
        title( display_out, 'Output Image' ) ;set( display_out, 'YTickLabel', [ ] ) ; set( display_out, 'XTickLabel', [ ] ) ;
        
        hold(axh,'off')

    end

    img_out = uint8(img_out);