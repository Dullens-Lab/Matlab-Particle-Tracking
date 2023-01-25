# Particle Tracking tailored for the Dullens Lab

If you want to contribute (recommended), ask to be added to the repo team and, 

`git clone https://github.com/Dullens-Lab/Matlab-Particle-Tracking`

If you only want to use the code, download the repo your preferred way.

## Tutorial

## Image Filtering, `bpass()`

Three step filtering starting with a boxcar filter to remove long scale variation. Boxcar filtered image, 'img_box' is then filtered with a gaussian filter, to remove pixel noise. Finally, 'img_gaus', has any pixel values below 'baseline' set to zero. Any step can be skipped with a 'false' argument.

`img_out = bpass( img, box_filter, gaus_filter, baseline, display )`

`img`:          'array' 2D array of image pixel values.

`box_filter`:   'bool' Set to 'true' for highpass filtering.
              Set to 'false' to skip boxcar filtering of the input image.

`gaus_filter`: 'int|bool' Characteristic length scale of noise in pixels.
              Set to 'true' to apply an appropriate gausian filter based on the input image.
              Or provide any positive integer for manual control of the gausian kernel. 
              Set to 'false' to skip gausian filtering of the input image.

`baseline`:     'int|bool' Reset any pixel values below 'baseline' to 0.
              Set to 'false' to skip. 
              An input of '0' will output the same result as 'false', but the code will scan the image for any values below 0.

returns:      'array' 2D array of filtered image pixel values.

Notes:        Performs a bandpass by convolving with an appropriate kernel. You can think of this as 
              a two part process. First, a lowpassed image is produced by convolving the original 
              with a gausian. Next, a second lowpassed image is produced by convolving the original 
              with a boxcar function. By subtracting the boxcar version from the gausian version, we
              are using the boxcar version to perform a highpass.


The Matlab Particle Tracking Code originally developed by Daniel Blair and Eric Dufresne, based on the IDL algorithms developed by David Grier, John Crocker and Eric Weeks.

Original Matlab code is distributed at https://site.physics.georgetown.edu/matlab/

Original IDL code is distributed at https://physics.emory.edu/faculty/weeks/idl/

