function SDBF = aux_fcn_SDBF(x, Up, Lp, k)
    % Smooth DeadBand Function with Numerical Stability
    %
    % x : input variable (scalar or vector)
    % Up: upper deadband limit
    % Lp: lower deadband limit
    % k : smoothing factor (positive scalar)

    % Midpoint and half-width
    alpha = (Up + Lp) / 2;
    mu    = (Up - Lp) / 2;
    eps_denom = 1e-12;

    % Normalized variable
    z = (x - alpha) ./ (mu + eps_denom);

    % Stable nonlinear term
    smooth = z ./ ((1 + abs(z).^k).^(1/k));

    % Final output
    SDBF = x - alpha - mu * smooth;

    % Clamp extreme results for safety
    SDBF = min(max(SDBF, Lp - alpha - mu), Up - alpha + mu);

    % Replace any undefined values with zero
    SDBF(isnan(SDBF) | isinf(SDBF)) = 0;
end
