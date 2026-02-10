function traj_zeroed = zero_trajectories( trajectories )
    
    x = trajectories(:,1);
    y = trajectories(:,2);
    t = trajectories(:,3);
    id = trajectories(:,4);
    
    % Zero all the trajectories so they start at the same place
    x0 = zeros(size(x));
    y0 = zeros(size(y));
    
    ids = unique(id);
    
    % Loop over each particle ID
    for k = 1:numel(ids)
        mask = (id == ids(k));
    
        % take the first time point of this trajectory
        idx = find(mask, 1, 'first');
    
        x_ref = x(idx);
        y_ref = y(idx);
    
        x0(mask) = x(mask) - x_ref;
        y0(mask) = y(mask) - y_ref;
    end
    
    traj_zeroed = [x0, y0, t, id];
end