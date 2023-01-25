# Particle Tracking tailored for the Dullens Lab

If you want to contribute (recommended), ask to be added to the repo team and, 

`git clone https://github.com/Dullens-Lab/Matlab-Particle-Tracking`

If you only want to use the code, download the repo your preferred way.

## Tutorial

### Image Filtering, `bpass()`

Three step filtering starting with a boxcar filter to remove long scale variation. Boxcar filtered image, 'img_box' is then filtered with a gaussian filter, to remove pixel noise. Finally, 'img_gaus', has any pixel values below 'baseline' set to zero. Any step can be skipped with a 'false' argument.

`img_out = bpass( img, box_filter, gaus_filter, baseline, display )`

`img` 2D array of image pixel values.

`box_filter` Set to `true` for highpass filtering. Set to `false` to skip boxcar filtering of the input image.

`gaus_filter` Characteristic length scale of noise in pixels. Set to `true` to apply an appropriate gaussian filter based on the input image. Or provide any positive integer for manual control of the gaussian kernel. Set to `false` to skip gaussian filtering of the input image. For either auto or manual values equal to 1 the image is assumed to be good enough to not need gaussian filtering and will be skipped.

`baseline` Reset any pixel values below 'baseline' to 0.Set to 'false' to skip. An input of '0' will output the same result as 'false', but the code will scan the image for any values below 0.

`img_out` 2D array of filtered image pixel values.

#### Ideal case example 

With an appropriate brightfield image, where the sample if focused to maximize the bright spot at it's centre but not saturating any pixels.

|![Ideal input image](/img/img_in_150.jpg)|

| Ideal brightfield image |

`img_in = imread( '../img/tutorial_150.tif' ) ;`
`[img_gaus, img_box, img_out] = bpass( img, true, 2, 120, true ) ;`

![Boxcar filtered image](/img/img_box_150.jpg)
| Boxcar filtered image |

![Boxcar then Gaussian filtered image](/img/img_gaus_150.jpg)
| Boxcar then Gaussian filtered image |

![Final output image](/img/img_out_150.jpg)
| Final output image |



The Matlab Particle Tracking Code originally developed by Daniel Blair and Eric Dufresne, based on the IDL algorithms developed by David Grier, John Crocker and Eric Weeks.

Original Matlab code is distributed at https://site.physics.georgetown.edu/matlab/

Original IDL code is distributed at https://physics.emory.edu/faculty/weeks/idl/

