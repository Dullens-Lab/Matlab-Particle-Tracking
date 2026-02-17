% trajectories: N x 4 array [x, y, t, id]
x  = traj_corr_um(:,2);
t  = traj_corr_um(:,3);
id = traj_corr_um(:,4);

dt = 1/fps;

% Sort by (id, t) just to be safe
[~, idx] = sortrows([id, t], [1 2]);
x  = x(idx);
id = id(idx);

% Find indices where a new particle starts
newParticle = [true; diff(id) ~= 0];

% Compute dx
dx_all = diff(x);

vel = dx_all / dt ;

% Remove diffs that jump between particles
dx_all(newParticle(2:end)) = [];

% --- Set bin size (you control this) ---
dx_bin = 20e-9;   % <-- set this to whatever makes sense in your units


% Histogram -> PDF
[p, edges] = histcounts(dx_all, 'BinWidth',dx_bin, 'Normalization', 'probability');
centers = 0.5*(edges(1:end-1) + edges(2:end));

% --- Plot ---
figure;
plot(centers, p, 'o');
xlabel('\Delta x');
ylabel('P(\Delta x)');
title('Step-size distribution over all particles');
grid on;
box on;

% Example data
x = centers(:);
y = p(:);

% Define Gaussian model: a*exp(-(x-b)^2/(2*c^2)) + d
ft = fittype('a*exp(-((x-b)^2)/(2*c^2)) + d', ...
             'independent', 'x', 'coefficients', {'a','b','c','d'});

% Initial guesses help convergence
start.a = max(y) - min(y);   % amplitude
start.b = x(y == max(y));    % mean
start.c = std(x);            % width (sigma)
start.d = min(y);            % offset

[fitresult, gof] = fit(x, y, ft, 'StartPoint', [start.a, start.b(1), start.c, start.d]);

% Extract parameters
A     = fitresult.a;
mu    = fitresult.b;
sigma = abs(fitresult.c);
offset= fitresult.d;

% Plot
figure;
plot(x, y, 'o'); hold on;
plot(fitresult, x, y);
legend('Data','Gaussian fit');
