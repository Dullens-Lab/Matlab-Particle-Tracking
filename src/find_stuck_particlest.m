function [stuckIDs, scores] = find_stuck_particlest(data, opts)
% data: (N x 4) [x y t id]
% Finds particles that do NOT follow ensemble drift (likely stuck).

arguments
    data (:,4) double
    opts.minPoints (1,1) double = 20
    opts.rMin (1,1) double = 50           % minimal net displacement to be considered "moving" (units of x,y)
    opts.sMin (1,1) double = 0.3          % drift-following strength threshold
end

x = data(:,1); y = data(:,2); t = data(:,3); id = data(:,4);

% --- COM trajectory per unique time ---
tu = unique(t);
xcom = zeros(size(tu));
ycom = zeros(size(tu));

for k = 1:numel(tu)
    idx = (t == tu(k));
    xcom(k) = mean(x(idx));
    ycom(k) = mean(y(idx));
end

% Map time -> com
% (use ismember for robust indexing)
[~, loc] = ismember(t, tu);
xcom_at = xcom(loc);
ycom_at = ycom(loc);

ids = unique(id);
stuck = false(size(ids));

scores = table('Size',[numel(ids) 6], ...
    'VariableTypes', {'double','double','double','double','double','double'}, ...
    'VariableNames', {'id','n','driftMag','netMag','s','p'});

% Global drift direction/magnitude over full movie
dCOM = [xcom(end)-xcom(1), ycom(end)-ycom(1)];
driftMag = norm(dCOM);

for j = 1:numel(ids)
    idx = (id == ids(j));
    n = sum(idx);

    scores.id(j) = ids(j);
    scores.n(j) = n;
    scores.driftMag(j) = driftMag;

    if n < opts.minPoints || driftMag < eps
        scores.netMag(j) = NaN;
        scores.s(j) = NaN;
        scores.p(j) = NaN;
        continue
    end

    % Use the particle's first/last time in its own track
    tj = t(idx);
    [t1, i1] = min(tj);
    [t2, i2] = max(tj);

    xj = x(idx); yj = y(idx);

    r1 = [xj(i1), yj(i1)];
    r2 = [xj(i2), yj(i2)];
    dR = r2 - r1;

    netMag = norm(dR);
    s = dot(dR, dCOM) / (driftMag^2);      % ~1 if follows drift, ~0 if stuck
    perp = norm(dR - s*dCOM) / driftMag;   % deviation from drift direction

    scores.netMag(j) = netMag;
    scores.s(j) = s;
    scores.p(j) = perp;

    % Classify stuck: doesn't move AND doesn't follow drift
    stuck(j) = (netMag < opts.rMin) && (s < opts.sMin);
end

stuckIDs = scores.id(stuck);
end