
function [trajectories_free, stuckIDs ] = find_stuck_particles(trajectories, minPoints, rMax)
    % trajectories: [ x y t id ]
    % minPoints:    minimum detections required per ID (default 100)
    % rMax:         max variance to be considered stuck [same units as x,y]
    
    if nargin < 2 || isempty( minPoints ), minPoints = 100; end
    if nargin < 3 || isempty( rMax ),      rMax      = 50; end
    
    % Seperate into x, y, t, id for easier handling
    x  = trajectories( :, 1 ) ; y  = trajectories( :, 2 ) ; t  = trajectories( :, 3 ) ; id = trajectories( :, 4 ) ;
    
    % Get the unique particle IDs
    ids = unique( id ) ;
    % Initialise a boolean array to keep track of which IDs are stuck
    stuck = false( size( ids ) ) ;
    
    % Loop over each unique ID and calculate the radius of gyration for its trajectory, x and y separately. If the minimum of the two is below the threshold, we consider the particle to be stuck.

    for j = 1:numel(ids)
        m = (id == ids(j));
        if nnz(m) < minPoints, continue; end
    
        tj = t(m);
        xj = x(m);
        yj = y(m);

        % Center of mass of the trajectory
        x_cm = mean( xj ) ;
        y_cm = mean( yj ) ;

        % Radius of Gyration in 1D for x and y separately
        Xg = sqrt( mean( ( xj - x_cm ).^2 ) ) ;
        Yg = sqrt( mean( ( yj - y_cm ).^2 ) ) ;

        stuck(j) = min( [Xg Yg] ) <= rMax;
    end
    
    stuckIDs = ids(stuck);

    trajectories_free = trajectories( ~ismember(trajectories(:,4), stuckIDs), : );

    
    figure;
    hold on;
    
    for k = 1:numel(stuckIDs)
        thisID = stuckIDs(k);
    
        mask = trajectories(:,4) == thisID;
    
        x = trajectories(mask,1);
        y = trajectories(mask,2);
        t = trajectories(mask,3);
    
        % sort by time in case rows are not ordered
        [~, idx] = sort(t);
        x = x(idx);
        y = y(idx);
    
        plot(x, y, 'o') ;
    end
    
    axis equal;
    xlabel('x (pixels)');
    ylabel('y (pixels)');
    xlim([0 2592]); ylim([0 1944])
    title('Trajectories of stuck particles');
    hold off;