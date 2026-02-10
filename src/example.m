%
% Particle Tracking Template
%

clearvars ; addpath(genpath('../')) % Add all "Matlab Particle Tracking" directories

% Load paramameters, make sure and edit them to suite your data
tutorial_params 

for frame = 1 : 1200
    
    % Load current image
    img_in = double( imread( [ '/Users/arrancurran/Documents/radboud/Teaching/Soft Matter Practical/Test Data/captured_images_2025-01-25_12-48-01/img_' num2str( frame ) '.tiff' ] ) ) ;
    
    % Filter the image
    filtered_image = bpass( img_in, false, 120, false ) ;
    
    % Find the peak pixels
    est_pks = pkfnd( filtered_image, 120, 31 ) ;
    
    % Calculate the sub pixel centroids
    cntrds = cntrd( filtered_image, est_pks, 31, true, frame ) ;
    
    % Collect the result
    centroids = [ centroids ; cntrds ] ;
end

% Pass x, y, frame number to track()
trajectories = track( centroids( :, [ 1 2 5 ] ), maxdisp, param ) ;

trajectories_um = trajectories ;
trajectories_um( :, [ 1 2 ] ) = trajectories_um( :, [ 1 2 ] ) / pxum ;

[ msd_2d, msd_x, msd_y, tau, msd_count ] = calcMSD( trajectories_um, fps ) ;

% Plot particle brightness vs size.
figure ; plot( centroids(:,3), centroids(:,4), 'o')
xlabel('Peak Brightness'); ylabel('Estimated Sizes');

% Plot X and Y
figure ; plot( trajectories(:,1), trajectories(:,2), '.')
xlabel('x (pixels)'); ylabel('y (pixels)');
xlim([0 2592]); ylim([0 1944])

% Plot 2D, and 1D MSDs
figure ; plot( tau, msd_2d, 'ko', tau, msd_x, 'ro', tau, msd_y, 'bo' )
xlabel('Lag time, \tau, (s)'); ylabel('\langle \Deltar( \tau )^2 \rangle (m^2)');

D_T = (k_B * T) / (6 * pi * eta * a ) ;

D_m = MSD_m / 2 ;


% Find stuck particles based on their trajectories
% [trajectories_free,stuckIDs] = find_stuck_particles(trajectories, [], 100);


% Removing drift
% traj_zeroed = zero_trajectories(trajectories_free);
% traj_mean = mean_trajectories(traj_zeroed) ;
% Plot mean x and mean y and fit lines to get m_x and m_y
% plot(traj_mean(:,3),traj_mean(:,1),'o')

% m_x = -0.8;
% m_y = 0.1055 ;


% traj_corr = remove_drift(trajectories_free, m_x,m_y);
% 
% traj_corr_um = traj_corr ;
% traj_corr_um( :, [ 1 2 ] ) = traj_corr_um( :, [ 1 2 ] ) / pxum ;
