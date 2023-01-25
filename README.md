# Particle Tracking tailored for the Dullens Lab

If you want to get updates from the repo or if you want to contribute (recommended), ask to be added to the repo team and, 

`git clone https://github.com/Dullens-Lab/Matlab-Particle-Tracking`

If you only want to use the code, download the repo your preferred way.

### Image Filtering `bpass()`

`img_out = bpass( img_in, hpass, lpass, baseline, display )`

Three step image manipulation starting with a high frequency pass filter to remove long scale variations. The high pass filtered image, `img_hpass` is then filtered with a low (gaussian) pass filter, to remove pixel noise. Finally, `img_lpass`, has any pixel values below `baseline` set to zero. Any step can be skipped with a `false` argument.

`img_in` 2D array of image pixel values.

`hpass` Set to `true` for highpass filtering. Set to `false` to skip.

`lpass` Set to `true` to apply a gaussian filter with a strength calculated from the input image. Provide any positive integer for manual control of the gaussian kernel. Set to `false` to skip. For either auto or manual, if the strength of the filter is equal to 1 the image is assumed to be good enough amd gaussian filtering and will be skipped.

`baseline` Resets any pixel values below `baseline` to 0. Set to 'false' to skip.

`display` Plot the image and pixel distribution at each stage of the filtering. Set to `false` or leave blank to skip.

Returns `img_out` 2D array of filtered image pixel values.

`img_hpass` and `img_lpass` can be returned with 

`[ img_out, img_hpass ] = bpass()` and 

`[ img_out, ~, img_lpass ] = bpass()`, respectively.


### Find Colloids Central Pixel `pkfnd()`

`est_pks = pkfnd( img, threshold, excl_dia )`

Find particle positions in an image with pixel accuracy. The output here is expected to be an argument in `cntrd()`. Rarely would the output here be sufficient for analysis of particle dynamics.

All non-zero pixels within the exclusion radius `excl_dia / 2` of the image edges are eliminated. Checks each pixel and the 8 nearest neighbors to see which is brightest. Each remaining pixel is then checked to see if it is the brightest in a region of interest defined by `excl_dia`. 

`img_in` 2D array of image pixel values. Ideally each colloid is represented by only a few non-zero pixels.
   
`threshold` Eliminates pixel values below `threshold`.
   
`excl_dia` Diameter, in pixels, over which to exclude all but the brightest pixel. Also excludes pixels within `excl_dia / 2` of image edges.

Returns `est_pks` N x 2 array containing, pixelated coordinates of local maxima.


### Calculate the Colloids Centroid `cntrd()`

`cntrds = cntrd( img_in, est_pks, excl_dia, apply_mask )`

Calculates the centroid of a colloids position to sub-pixel accuracy.

`img_in` 2D array of image pixel values. Ideally each colloid is represented by as many pixels as possible.

`est_pks` Coordinates of central pixel for each of colloids.

`excl_dia` Diameter, in pixels, over which to calculate the centroid.

Returns `cntrds` N x 4 array containing, each colloids centroid, the brightest pixel and estimated radius.

    `cntrds(:,1)` x-coordinates
    `cntrds(:,2)` y-coordinates
    `cntrds(:,3)` brightest pixel
    `cntrds(:,4)` estimated radii

Note that sub-pixel accuracy is dependent on the number of pixels over which the centroid is calculated. To check for pixel bias, plot a histogram of the fractional parts of the resulting locations. For example,

`hist( cntrds( :, 1 ) ./ floor( cntrds( :, 1 ) ) )`.



#### Ideal case example 

With an appropriate brightfield image, where the sample if focused to maximize the bright spot at it's centre but not saturating any pixels.

|![Ideal input image](/img/img_in_150.jpg)|
|:--:|
| Ideal brightfield image |

`img_in = imread( '../img/tutorial_150.tif' ) ;`

`[img_lpass, img_hpass, img_out] = bpass( img, true, 2, 120, true ) ;`

|![Boxcar filtered image](/img/img_hpass_150.jpg)|
|:--:|
| `img_hpass` |

The boxcar filter is convolution with a 3x3 array 

|![Boxcar then gaussian filtered image](/img/img_lpass_150.jpg)|
|:--:|
| `img_lpass` |

|![Final output image](/img/img_out_150.jpg)|
|:--:|
| `img_out`|



The Matlab Particle Tracking Code originally developed by Daniel Blair and Eric Dufresne, based on the IDL algorithms developed by David Grier, John Crocker and Eric Weeks.

Original Matlab code is distributed at https://site.physics.georgetown.edu/matlab/

Original IDL code is distributed at https://physics.emory.edu/faculty/weeks/idl/

