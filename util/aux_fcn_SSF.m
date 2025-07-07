function SSF = aux_fcn_SSF(x, U, L, k)
    % Smooth Saturation Function with Numerical Stability
    %
    % x: input state (scalar or vector)
    % U: upper bound
    % L: lower bound
    % k: smoothness factor (positive scalar)

    % Parameters
    alpha = (U + L) / 2;
    mu    = (U - L) / 2;
    eps_denom = 1e-12;

    % Core expression
    z = (x - alpha) ./ (mu + eps_denom);        % normalized input
    smooth = z ./ ((1 + abs(z).^k).^(1/k));     % stabilized nonlinear function

    % Final SSF value
    SSF = alpha + mu * smooth;

    % Optional: clamp output to bounds for safety
    SSF = min(max(SSF, L), U);

    % Catch any remaining numerical problems
    SSF(isnan(SSF) | isinf(SSF)) = alpha;  % reset to midpoint if undefined
end
