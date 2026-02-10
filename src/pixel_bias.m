function out = pixel_bias(trajectories, nbins, doPlot)
%PIXEL_BIAS_HIST  Simple histogram-based pixel bias check.
%
%   out = pixel_bias_hist(x, y)
%   out = pixel_bias_hist(x, y, nbins, doPlot)
%
% Inputs:
%   trajectories   : vectors of tracked positions (in pixels)
%   nbins  : number of bins in [0,1) (default: 20)
%   doPlot : true/false to plot histograms (default: true)
%
% Output struct:
%   out.biasX, out.biasY : L1 distance from uniform (0 = perfectly uniform)
%   out.histX, out.histY : histogram counts
%   out.fx, out.fy       : fractional coordinates
%
% Interpretation:
%   - If biasX or biasY ~ 0 → no obvious pixel bias
%   - Larger values → stronger pixel locking / bias
%   - Look at the plots: peaks at specific fractions (0, 0.5, etc.) = trouble

    if nargin < 3 || isempty(nbins)
        nbins = 40;
    end
    if nargin < 4 || isempty(doPlot)
        doPlot = true;
    end
    
    x = trajectories(:,1);
    y = trajectories(:,2);

    % Keep only finite points
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    % Fractional parts in [0,1)
    fx = mod(x, 1);
    fy = mod(y, 1);

    edges = linspace(0, 1, nbins+1);

    histX = histcounts(fx, edges);
    histY = histcounts(fy, edges);

    % Normalize to probabilities
    px = histX / sum(histX);
    py = histY / sum(histY);

    % Uniform distribution
    u = ones(1, nbins) / nbins;

    % Simple bias metric: L1 distance from uniform
    biasX = 0.5 * sum(abs(px - u));
    biasY = 0.5 * sum(abs(py - u));

    % Pack output
    out = struct();
    out.biasX = biasX;
    out.biasY = biasY;
    out.histX = histX;
    out.histY = histY;
    out.fx = fx;
    out.fy = fy;

    if doPlot
        centers = (edges(1:end-1) + edges(2:end)) / 2;

        figure;
        subplot(1,2,1);
        bar(centers, px, 1);
        xlim([0 1]);
        xlabel('Fractionals in x'); ylabel('Probability Distribution');
        title(sprintf('Estimated Bias = %.3f', biasX));

        subplot(1,2,2);
        bar(centers, py, 1);
        xlim([0 1]);
        xlabel('Fractionals in y'); ylabel('Probability Distribution');
        title(sprintf('Estimated Bias = %.3f', biasY));
    end
end
