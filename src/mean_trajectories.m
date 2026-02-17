function traj_mean = mean_trajectories( traj_zeroed )
    
    x = traj_zeroed(:,1);
    y = traj_zeroed(:,2);
    t = traj_zeroed(:,3);

    % Average over particles to get one trajectory
    tvals = unique(t);
    
    x_mean = zeros(size(tvals));
    y_mean = zeros(size(tvals));
    
    % Loop over all frames
    for k = 1:numel(tvals)
        % Get all indicies with the same t
        mask = (t == tvals(k));
        x_mean(k) = mean(x(mask), 'omitnan');
        y_mean(k) = mean(y(mask), 'omitnan');
    end
    
    traj_mean = [x_mean, y_mean, tvals];
end