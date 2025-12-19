% 
% Calculate the Mean Squared Displacement (MSD) of multiple tracked particles over time 
% and plot the result.
% 
% MSD( dt ) = < [ x( t + dt ) - x( t ) ]^2 >
% 
%   `msd = calcMSD( Particle )`
%   
%   `Particle`: A matrix containing particle tracking data with four columns:
%       
%       `x` — x-coordinate of the particle’s position.
%       `y` — y-coordinate of the particle’s position.
%       `t` — Time/frame number of the observation.
%       `ID` — A unique identifier for each particle.
%       
%   `msd`: A two-column matrix where:
%   
%   `tau` — Time lag between observations.
%   `MSD` — Mean squared displacement for the corresponding time lag.
% 
% Overview:
% 
% 1. Determine the maximum recorded time (`lmax`) and list of unique particle IDs.
% 
% 2. Initialize two arrays, `MSD` and `MSD_count`, to accumulate displacement sums and counts for each time lag.
% 
% 3. Loop over each particle ID:
%    - Extract the data for the current particle.
%    - For each possible time difference (`dt`), loop over all possible starting points within the particle’s trajectory.
%    - Calculate the real time difference (`realdt`) between the two frames, accounting for possible missing frames.
%    - Compute the squared displacement between the two positions (`v1` and `v2`).
%    - Accumulate the squared displacement in the `MSD` array and increment the count in `MSD_count`.
% 
% 4. Normalize the `MSD` array by the count of observations.
% 
% 5. Plot the MSD against time lag (`tau`).

function [msd, tau, msd_count] = calcMSD( Particle )
    
    % Length of trajectories
    lmax = max(Particle(:, 3));
    
    % Array of all particle IDs
    IDs = unique(Particle(:, 4));
    
    msd = zeros(lmax, 1);
    msd_count = zeros(lmax, 1);
    
    % Loop over each particle
    for i = 1:length(IDs)
        
        % Current particle ID
        thisID = IDs(i);
        
        % Current particle
        thisParticle = Particle(Particle(:, 4) == thisID, :);
        
        % Loop over current particle
        for dt = 0:size(thisParticle, 1) - 1
            
            valid_indices = (1:size(thisParticle, 1) - dt)';
            
            t1 = valid_indices;
            t2 = t1 + dt;
            
            % Calculate real dt (in case of missing frame)
            realdt = thisParticle(t2, 3) - thisParticle(t1, 3);
            
            % Calculate squared displacement
            displacement = sum((thisParticle(t2, 1:2) - thisParticle(t1, 1:2)).^2, 2);
            
            % Accumulate results
            for j = 1:length(realdt)
                msd_count(realdt(j) + 1) = msd_count(realdt(j) + 1) + 1;
                msd(realdt(j) + 1) = msd(realdt(j) + 1) + displacement(j);
            end
        end
    end
    
    msd = msd ./ msd_count;
    tau = 1:lmax;

end


