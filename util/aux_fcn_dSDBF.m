function dSDBF = aux_fcn_dSDBF(x, Up, Lp, k)
    % Derivative of the Smooth Deadband Function
    % SDBF(x) = x - alpha - mu * ( z / (1 + z^k)^(1/k) )
    % dSDBF = 1 - d/dx[ z / (1 + z^k)^(1/k) ]

    eps_mu = 1e-6;  % to avoid division by zero
    alpha = 0.5 * (Up + Lp);
    mu    = max(0.5 * (Up - Lp), eps_mu);

    % Normalize
    z = (x - alpha) / mu;
    z = max(min(z, 100), -100);  % clamp for numerical safety

    z_pow_k = z.^k;
    denom = (1 + z_pow_k);
    denom_root_k = denom.^(1/k);
    denom_root_kp1 = denom.^(1 + 1/k);

    % Compute dphi/dz
    dphi_dz = (1 ./ denom_root_k) - ((z .* z_pow_k) ./ denom_root_kp1);

    % Final derivative
    dSDBF = 1 - dphi_dz;

    % Safety cleanup
    dSDBF(~isfinite(dSDBF)) = 0;
end
