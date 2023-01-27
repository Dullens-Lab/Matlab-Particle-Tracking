# Particle Tracking tailored for the Dullens Lab

Find positions within a microscope image or series of images of single or multiple colloids.

If you want to get updates from this repo or if you want to contribute (recommended), ask to be added to the repo team and, 

`git clone https://github.com/Dullens-Lab/Matlab-Particle-Tracking`

If you only want to use the code, download the repo your preferred way.

### `bpass()` Image Filtering

`img_out = bpass( img_in, hpass, lpass, bckgrnd, display )`

Three step image manipulation starting with a high pass filter to remove long scale variations. The high pass filtered image, `img_hpass` is then filtered with a low pass filter, to remove pixel noise. Finally, `img_lpass`, has any pixel values below `backgrnd` set to zero. Any step can be skipped with a `false` argument.

`img_in` 2D array of image pixel values.

`hpass` Set to `true` for high pass filtering. Set to `false` to skip.

`lpass` Set to `true` to apply a gaussian filter with a strength calculated from `img_in`. Provide any positive integer for manual control of the gaussian kernel. Set to `false` to skip. For either `true` or `int value`, if the strength of the filter is equal to 1 the image is assumed to be good enough and `lpass` and will be skipped.

`backgrnd` Resets any pixel values below `backgrnd` to 0. Set to 'false' to skip.

`display` Display `img_in`, `img_hpass`, `img_lpass` and `img_out`. And calculates pixel distribution at each stage of the filtering. This can be usefully when optimising your input arguments. Set to `false` or leave blank to skip.

Returns `img_out` 2D array of filtered image pixel values.

`img_hpass` and `img_lpass` can be returned with 

`[ img_out, img_hpass ] = bpass()` and 

`[ img_out, ~, img_lpass ] = bpass()`, respectively.


### `pkfnd()` Find Colloids Central Pixel

`est_pks = pkfnd( img, threshold, excl_dia )`

Find colloid positions in an image with pixel accuracy. The output here is expected to be an argument in `cntrd()`. Rarely would the output here be sufficient for analysis of particle dynamics.

All non-zero pixels within the exclusion radius `floor( excl_dia / 2 )` of the image edges are eliminated. Checks each pixel and eliminates all but the brightest within a 3x3 array centered on the current pixel. After which remaining pixels are then checked to see if it is the brightest in a n x n array centered on the current pixel where `n = excl_dia`.

`img_in` 2D array of image pixel values. Ideally each colloid is represented by only a few non-zero pixels.
   
`threshold` Eliminates pixel values below `threshold`.
   
`excl_dia` Diameter, in pixels, over which to exclude all but the brightest pixel. Also excludes pixels within `excl_dia / 2` of image edges.

Returns `est_pks` N x 2 array containing, pixelated coordinates of local maxima.


### `cntrd()` Calculate the Colloid Centroids

`cntrds = cntrd( img_in, est_pks, excl_dia, apply_mask )`

Calculates the centroid of a colloid to sub-pixel accuracy.

`img_in` 2D array of image pixel values. Ideally each colloid is represented by as many pixels as possible.

`est_pks` Coordinates of central pixel for each colloid.

`excl_dia` Diameter, in pixels, over which to calculate the centroid.

Returns `cntrds` N x 4 array containing, centroids, the corresponding brightest pixel and estimated radius.

`cntrds(:,1)` x-coordinates

`cntrds(:,2)` y-coordinates

`cntrds(:,3)` brightest pixel

`cntrds(:,4)` estimated radii

Note that sub-pixel accuracy is dependent on the number of pixels over which the centroid is calculated. To check for pixel bias, plot a histogram of the fractional parts of the resulting locations. For example,

`hist( cntrds( :, 1 ) ./ floor( cntrds( :, 1 ) ) )`.


## Examples

### Ideal case

With an appropriate brightfield image, where monodisperse colloids are focused to maximize the contrast between colloids and the background but not saturating any pixels.

|![Ideal input image](/img/img_in_ideal.jpg)|  
|:--:|
| `img_in` |

In the above image there is negligible noise and negligible variation across the dimensions of the image. This means we can skip the long pass and high pass filtering steps and just apply an appropriate `backgrnd` to remove all of the background whilst retaining as many as possible pixels per colloid.

`img_out = bpass( img_in, false, false, 65 ) ;`

|![Ideal Output Image](/img/img_out_ideal.jpg)|
|:--:|
| `img_out` |

Now we can estimate the coordinates of each colloid to pixel level accuracy with `pkfnd()`.

`est_pks = pkfnd( img_out, 172, 11 )`

With `pkfnd()` we can apply a much higher threshold than with `bpass()` so as to reduce the number of candidate pixels per colloid. You can also check this by returning the total number of candidate pixels after the applied threshold with,

`[ est_pks, input_pk_pxs] = pkfnd( img_out, 172, 9 )`

Ideally `input_pk_pxs` is minimized without loosing any colloids.

Finally, we pass the `est_pks` to `cntrd()` along with either our filtered image `img_out`, or in this case, since we have such a good image we can simple pass the original image and apply the circular mask function,

`cntrds = cntrd( img_in, est_pks, 11, true )`

|![Ideal Output Image](/img/img_result_ideal.jpg)|
|:--:|
| `img_result` |

The above image shows the original image with the centroid markers (green crosshair), the `excl_dia` (green square) and the estimated radii (red circles).

It would normally be advisable to check for pixel biasing but since we only have a handful of particles we can simple just look at the `cntrds` values.


### Dealing with noise

Consider the case where we have a poor camera or low light levels (fluorescence for example). Images here will have high frequency noise that will get in the way of accurately calculating the centroids.

|<img src="/img/img_in_noisy.jpg" width="150">|
|:--:|
| `img_in` |



The Matlab Particle Tracking Code originally developed by Daniel Blair and Eric Dufresne, based on the IDL algorithms developed by David Grier, John Crocker and Eric Weeks.

Original Matlab code is distributed at https://site.physics.georgetown.edu/matlab/

Original IDL code is distributed at https://physics.emory.edu/faculty/weeks/idl/

