k_B         = 1.380649e-23; % J K^{-1}
T           = 298 ;
eta         = 0.89e-3 ;

centroids   = [] ;

a           = 1e-6 ;

fps         = 2 ;   % Image Acquistion frame rate

pxum        = 8 ;   % Microscope calibartion, pixels per micrometre [ /um ]
pxum        = pxum / 1e-6 ; % pixels per metre [ /m ]

excl_dia    = 31 ;  % Diameter, in pixels, where only one peak pixel will be recorded
backgrnd    = 120 ; % Background threshold
maxdisp     = 21 ;  % Estimated maximum displacement a colloid may undergoe between frames

param.mem   = 2 ;  % Number of frames to keep a lost particle in the memoroy
param.good  = 600 ; % At the end of track(), colloids with centroids < param.good are discarded 
param.dim   = 2 ;
param.quiet = 0 ;

excl_rad    = floor( excl_dia / 2 ) ;