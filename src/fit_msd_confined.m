% Simple confined MSD fit:
%   MSD(tau) = A * (1 - exp(-tau / tau_c))
%
% Usage:
%   [A, tau_c, msd_fit] = fit_msd_confined(tau, msd)
%
function [A, tau_c, msd_fit] = fit_msd_confined(tau, msd)

    tau = tau(:);
    msd = msd(:);

    valid = isfinite(tau) & isfinite(msd) & (tau >= 0);
    tau = tau(valid);
    msd = msd(valid);

    fit_fun = @(p, t) p(1) .* (1 - exp(-t ./ p(2)));

    A0 = max(msd);
    tau_pos = tau(tau > 0);
    if isempty(tau_pos)
        tau_c0 = 1;
    else
        tau_c0 = median(tau_pos);
    end

    p0 = [A0, tau_c0];
    obj_fun = @(p) sum((msd - fit_fun(abs(p), tau)).^2);
    p = abs(fminsearch(obj_fun, p0, optimset('Display', 'off')));

    A = p(1);
    tau_c = p(2);
    msd_fit = fit_fun(p, tau);

    figure;
    plot(tau, msd, 'ko', tau, msd_fit, 'r-', 'LineWidth', 1.5);
    xlabel('Lag time, tau (s)');
    ylabel('MSD');
    legend('Data', 'Fit', 'Location', 'best');
    grid on;
end
