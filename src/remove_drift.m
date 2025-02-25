function [correctedCOORD, comTraj] = remove_drift(COORD)
%REMOVE_DRIFT Removes drift from a structured particle tracking dataset.
%
%   INPUT:
%       COORD - 1D structure array where COORD(k) contains an (N x 6) array:
%               Columns: [X, Y, brightnes, radius, ?, Time], where X and Y are in cols 1 & 2, Time is in col 6.
%
%   OUTPUT:
%       correctedCOORD - Same structure as COORD but with drift-corrected X, Y positions.
%       comTraj       - (T x 3) array containing [Time, X_com, Y_com] for each frame.
%
%   USAGE:
%       correctedCOORD = remove_drift_struct(COORD);
%

    T = length(COORD);  % Number of frames

    % Initialize arrays for storing COM trajectory
    comTraj = zeros(T, 3);  % [Time, X_com, Y_com]
    
    % Compute center-of-mass (COM) trajectory
    for k = 1:T
        frameData = COORD(k).coordinates; % N x 6 array
        x_com = mean(frameData(:,1)); % Mean X
        y_com = mean(frameData(:,2)); % Mean Y
        time_k = frameData(1,6); % Time from first particle (assuming consistent time)
        
        comTraj(k, :) = [time_k, x_com, y_com]; % Store COM trajectory
    end

    % Subtract COM motion from each frame
    correctedCOORD = COORD; % Copy structure to modify
    for k = 1:T
        correctedCOORD(k).coordinates(:,1) = COORD(k).coordinates(:,1) - comTraj(k,2); % X corrected
        correctedCOORD(k).coordinates(:,2) = COORD(k).coordinates(:,2) - comTraj(k,3); % Y corrected
    end
end