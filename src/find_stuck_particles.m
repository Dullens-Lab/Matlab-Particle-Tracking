function stuckIDs = find_stuck_particles(trajectories, minPoints, rMax)
% stuckIDs = find_stuck_particles(trajectories, minPoints, rMax)
%
% trajectories: (N x 4) [x y t id]
% Returns IDs whose net displacement (first->last) is <= rMax.
%
% minPoints: minimum detections required per ID (default 20)
% rMax:      max net displacement to be considered stuck (default 50) [same units as x,y]

if nargin < 2 || isempty(minPoints), minPoints = 20; end
if nargin < 3 || isempty(rMax),      rMax      = 50; end

x  = trajectories(:,1);
y  = trajectories(:,2);
t  = trajectories(:,3);
id = trajectories(:,4);

ids = unique(id);
stuck = false(size(ids));

for j = 1:numel(ids)
    m = (id == ids(j));
    if nnz(m) < minPoints, continue; end

    tj = t(m);
    xj = x(m);
    yj = y(m);

    [~, i1] = min(tj);
    [~, i2] = max(tj);

    dR = [xj(i2)-xj(i1), yj(i2)-yj(i1)];
    stuck(j) = (hypot(dR(1), dR(2)) <= rMax);
end

stuckIDs = ids(stuck);
end

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



end