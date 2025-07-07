function dSDBF_dbd2 = aux_fcn_dSDBF_dbd2(x, U, L, k)
    % Derivative of the smooth deadband function with respect to U (dbd2)
    
    % Safety threshold to prevent division by near-zero
    eps_mu = 1e-6;

    % Compute alpha and mu safely
    alpha = 0.5 * (U + L);
    mu = max(0.5 * (U - L), eps_mu);  % prevent mu ~ 0

    % Safe z
    z = (x - alpha) / mu;

    % Compute stable phi and dphi_dz
    z_pow_k = z^k;
    denom = 1 + z_pow_k;
    denom_root_k = denom^(1/k);
    denom_root_kp1 = denom^(1 + 1/k);

    phi = z / denom_root_k;

    % Avoid div-by-zero in derivative
    dphi_dz = (1 + z_pow_k)^(-1/k);
    if abs(denom_root_kp1) > eps_mu
        dphi_dz = dphi_dz - (z_pow_k * z) / denom_root_kp1;
    end

    % Derivatives of alpha and mu w.r.t. U
    dalpha_dU = 0.5;
    dmu_dU = 0.5;

    % Derivatives of z
    dz_dalpha = -1 / mu;
    dz_dmu = -(x - alpha) / (mu^2);

    % Derivatives of phi
    dphi_dalpha = dphi_dz * dz_dalpha;
    dphi_dmu = dphi_dz * dz_dmu;

    % Final expression with safeguards
    dSDBF_dbd2 = -dalpha_dU ...
                 - dmu_dU * phi ...
                 - mu * (dphi_dalpha * dalpha_dU + dphi_dmu * dmu_dU);

    % Clip any numerical noise
    if ~isfinite(dSDBF_dbd2)
        dSDBF_dbd2 = 0;
    end
end

