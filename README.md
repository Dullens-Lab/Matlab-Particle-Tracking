# Particle Tracking tailored for the Dullens Lab

If you want to contribute (recommended), ask to be added to the repo team and, 

`git clone https://github.com/Dullens-Lab/Matlab-Particle-Tracking`

If you only want to use the code, download the repo your preferred way.

## Tutorial

### Image Filtering, `bpass()`

Three step image manipulation starting with a high frequency pass filter to remove long scale variations. The high pass filtered image, `img_hpass` is then filtered with a low (gaussian) pass filter, to remove pixel noise. Finally, `img_lpass`, has any pixel values below `baseline` set to zero. Any step can be skipped with a `false` argument.

`img_out = bpass( img_in, hpass, lpass, baseline, display )`

`img_in` 2D array of image pixel values.

`hpass` Set to `true` for highpass filtering. Set to `false` to skip.

`lpass` Set to `true` to apply a gaussian filter with a strength calculated from the input image. Provide any positive integer for manual control of the gaussian kernel. Set to `false` to skip. For either auto or manual, if the strength of the filter is equal to 1 the image is assumed to be good enough amd gaussian filtering and will be skipped.

`baseline` Reset any pixel values below `baseline` to 0. Set to 'false' to skip.

`display` Plot the image and pixel distribution at each stage of the filtering. Set to `false` or leave blank to skip.

`img_out` 2D array of filtered image pixel values.

`img_hpass` and `img_lpass` can be returned with `[ img_out, img_hpass ] = bpass()` and `[ img_out, ~, img_lpass ] = bpass()`, respectively.


### Find peak pixels, `pkfnd()`







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

