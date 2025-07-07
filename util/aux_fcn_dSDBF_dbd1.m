function dSDBF_dbd1 = aux_fcn_dSDBF_dbd1(x, U, L, k)
    % Derivative of the smooth deadband function with respect to L (dbd1)

    % Small epsilon to avoid division by near-zero mu
    eps_mu = 1e-6;

    % Safely compute alpha and mu
    alpha = 0.5 * (U + L);
    mu = max(0.5 * (U - L), eps_mu);  % avoid mu ≈ 0

    % Normalized variable
    z = (x - alpha) / mu;

    % Compute stable phi and its derivative
    z_pow_k = z^k;
    denom = 1 + z_pow_k;
    denom_root_k = denom^(1/k);
    denom_root_kp1 = denom^(1 + 1/k);

    phi = z / denom_root_k;

    % Derivative of phi w.r.t. z
    dphi_dz = (1 + z_pow_k)^(-1/k);
    if abs(denom_root_kp1) > eps_mu
        dphi_dz = dphi_dz - (z_pow_k * z) / denom_root_kp1;
    end

    % Derivatives of alpha and mu w.r.t. L
    dalpha_dL = 0.5;
    dmu_dL = -0.5;

    % Chain rule for z
    dz_dalpha = -1 / mu;
    dz_dmu = -(x - alpha) / mu^2;

    % Chain rule for phi
    dphi_dalpha = dphi_dz * dz_dalpha;
    dphi_dmu = dphi_dz * dz_dmu;

    % Full derivative expression
    dSDBF_dbd1 = -dalpha_dL ...
                 - dmu_dL * phi ...
                 - mu * (dphi_dalpha * dalpha_dL + dphi_dmu * dmu_dL);

    % Final sanity clip
    if ~isfinite(dSDBF_dbd1)
        dSDBF_dbd1 = 0;
    end
end
