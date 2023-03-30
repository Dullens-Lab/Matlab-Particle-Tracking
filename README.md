# Particle Tracking tailored for the Dullens Lab

Find positions within a microscope image or series of images of single or multiple colloids.

If you want to get updates from this repo or if you want to contribute (recommended), ask to be added to the repo team and, 

`git clone https://github.com/Dullens-Lab/Matlab-Particle-Tracking`

## Jump To

- [Functions](#functions)
    * [Image Preparation with bpass()](#image-filtering)
    * [Finding Peak Pixels with pkfnd()](#find-peak-pixel)
    * [Sub-pixel Centroids with cntrd()](#calculate-the-colloid-centroids)

- [Usage](#usage)
    * [Ideal Case](#ideal-case)
    * [Pixel Biasing](#pixel-biasing)

## Functions

### Image Filtering

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


### Find Peak Pixel

`est_pks = pkfnd( img, threshold, excl_dia )`

Find colloid positions in an image with pixel accuracy. The output here is expected to be an argument in `cntrd()`. Rarely would the output here be sufficient for analysis of particle dynamics.

All non-zero pixels within the exclusion radius `floor( excl_dia / 2 )` of the image edges are eliminated. Checks each pixel and eliminates all but the brightest within a 3x3 array centered on the current pixel. After which remaining pixels are then checked to see if it is the brightest in a n x n array centered on the current pixel where `n = excl_dia`.

`img_in` 2D array of image pixel values. Ideally each colloid is represented by only a few non-zero pixels.

`threshold` Eliminates pixel values below `threshold`.

`excl_dia` Diameter, in pixels, over which to exclude all but the brightest pixel. Also excludes pixels within `excl_dia / 2` of image edges.

Returns `est_pks` N x 2 array containing, pixelated coordinates of local maxima.

### Calculate the Colloid Centroids

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


## Usage

### Ideal Case

A brightfield image, where monodisperse colloids are focused to maximize the contrast between colloids centre and the background, but not saturating any pixels.

|![Ideal input image](/img/img_in_ideal.jpg)|  
|:--:|
| `img_in` |

In the above image there is negligible noise and negligible variation across the image. We can skip the long pass and high pass filtering steps and just apply a `backgrnd` to zero all of the background pixels whilst retaining as many as possible pixels per colloid.

`img_out = bpass( img_in, false, false, 65 ) ;`

|![Ideal Output Image](/img/img_out_ideal.jpg)|
|:--:|
| `img_out` |

Now we can find the coordinates of the brightest pixel for each colloid with `pkfnd()`.

`est_pks = pkfnd( img_out, 172, 11 )`

With `pkfnd()` we can apply a much higher threshold than with `bpass()` so as to reduce the number of candidate pixels per colloid. You can also check this by returning the total number of candidate pixels after the applied threshold with,

`[ est_pks, input_pk_pxs] = pkfnd( img_out, 172, 9 )`

Ideally `input_pk_pxs` is minimized without loosing any colloids.

Finally, we pass the `est_pks` to `cntrd()` along with our filtered image `img_out`, an exclusion diameter in pixels, and a boolean,

`cntrds = cntrd( img_in, est_pks, 11, true )`

|![Ideal Output Image](/img/img_result_ideal.jpg)|
|:--:|
| `img_result` |

The above image shows the original image with the centroid markers (green crosshair), the `excl_dia` (green square) and the estimated radii (red circles).

### Pixel Biasing

It is advisable to check for pixel biasing by looking at the fractional components of the resulting centroids,

`hist( mod( cntrd( :, 1 ), 1 ) )` for example.

For a significant number of centroids, the distribution of the fractional component of the centroids should have a uniform distribution.


### Dealing with noise

If you have a particularly noisy image, the high frequency noise will propagate through to the calculated centroid.

Try running `noise_on_pos.m` found in `test/` where an image of a colloid with a random, uniform distribution of pixel noise is used to calculate the centroid with and without filtering the image with a low pass filter.

An example image with noise and the resulting filtered image with `bpass( img_noise, false, 5, 80)` is shown below.

|![Noisy Image](/img/img_in_noisy.jpg) ![Filtered Image](/img/img_in_noisy_filtered.jpg)|
|:--:|
| Noisy image and filtered image. |

To demonstrates the effect of noise, the distribution of the calculated centroids for 10<sup>5</sup> images with random uniform noise is shown below.

|![Centroid distribution](/img/cntrd_dist.jpg)|
|:--:|
| Centroid distribution for unfiltered and filtered images. |

#### Sub-Pixel Biasing

Aggressive filtering with the low pass filter can result in sub-pixel biasing due to [ringing artifacts](https://en.wikipedia.org/wiki/Ringing_artifacts). The larger the `lpass` argument the higher the suppression of the high frequencies variations and the image becomes bandlimited, leading to ringing at points of high contrast (ie from background pixels to colloid pixels).

|![Ringing Artifacts](/img/ringing_artifacts.jpg=100x100 ) ![Resulting Centroids Distribution](/img/ringing_artifacts_hist.jpg)|
|:--:|
| Centroids resulting from ringing artifacts and resulting distribution. |

### Credits

The Matlab Particle Tracking Code was originally developed by Daniel Blair and Eric Dufresne, based on the IDL algorithms developed by David Grier, John Crocker and Eric Weeks.

Original Matlab code is distributed at https://site.physics.georgetown.edu/matlab/

Original IDL code is distributed at https://physics.emory.edu/faculty/weeks/idl/