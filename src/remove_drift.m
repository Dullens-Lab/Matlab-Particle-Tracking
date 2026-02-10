function traj_corrected = remove_drift( trajectories, m_x, m_y)
    
    x = trajectories(:,1);
    y = trajectories(:,2);
    t = trajectories(:,3);
    id = trajectories(:,4);

    % After fitting to meanTraj remove drift from trajectories
    
    % Drift slope (units: xy-units per t-unit)
    
    
    % Compute drift term
    x_drift = m_x * t;
    y_drift = m_y * t;
    
    % Subtract drift
    x_corr = x - x_drift;
    y_corr = y - y_drift;
    
    % Reassemble corrected data
    traj_corrected = [x_corr, y_corr, t, id];

end