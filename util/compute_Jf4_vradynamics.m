function Jf4 = compute_Jf4_vradynamics(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    Jf4 = zeros(10, 10);
    eps_denom = 1e-6;
    idx = max(k-1, 1);

    % --- RK4 state update ---
    xk4 = xk;
    xk4(1:5) = xk(1:5) + k3(1:5);

    % --- Extract safe parameters ---
    Trv   = max(xk(6), eps_denom);
    kqv   = max(xk(7), eps_denom);
    Tg    = max(xk(8), eps_denom);
    Tiq   = max(xk(9), eps_denom);
    Tpord = max(xk(10), eps_denom);

    % --- Inputs ---
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    Vf    = max(xk4(1), eps_denom);

    % === Row 1: d(x1)/dt ===
    Jf4(1,1) = -dt / Trv;
    Jf4(1,6) = dt * (xk4(1) - Vt) / Trv^2;

    % === Row 2: d(x2)/dt ===
    Jf4(2,1) = dt * Qref / (Tiq * Vf^2);
    Jf4(2,2) = -dt / Tiq;
    Jf4(2,9) = dt * (xk4(2) - Qref / Vf) / Tiq^2;

    % === Row 3: d(x3)/dt — Keep kqv sensitivity, fixed dbds ===
    deltaV     = p.Vref - xk4(1);
    SDBF       = aux_fcn_SDBF(deltaV, p.dbd2, p.dbd1, p.k);
    dSDBF_dV   = -aux_fcn_dSDBF(deltaV, p.dbd2, p.dbd1, p.k);

    Iqv        = kqv * SDBF;
    dIqv_dVt   = kqv * dSDBF_dV;
    dIqv_dkqv  = SDBF;

    SSF_Iqv    = aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k);
    dSSF_Iqv   = aux_fcn_dSSF(Iqv, p.Iqh1, p.Iql1, p.k);

    psi        = xk4(2) - SSF_Iqv;
    dSSF_psi   = aux_fcn_dSSF(psi, p.Iqmax, p.Iqmin, p.k);

    Jf4(3,1) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dVt) / Tg;
    Jf4(3,2) = dt * dSSF_psi / Tg;
    Jf4(3,3) = -dt / Tg;
    Jf4(3,7) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dkqv) / Tg;
    Jf4(3,8) = dt * (xk4(3) - aux_fcn_SSF(psi, p.Iqmax, p.Iqmin, p.k)) / Tg^2;

    % === Row 4: d(x4)/dt ===
    SSF_P   = aux_fcn_SSF(xk4(4), p.Pmax, p.Pmin, p.k);
    dSSF_P  = aux_fcn_dSSF(xk4(4), p.Pmax, p.Pmin, p.k);
    Jf4(4,4)   = dt * (-1/Tpord * dSSF_P - p.kw * (dSSF_P - 1)) - dt / Tpord;
    Jf4(4,10)  = dt * (xk4(4) - Pref) / Tpord^2;

    % === Row 5: d(x5)/dt ===
    SSF_Px4 = aux_fcn_SSF(xk4(4), p.Pmax, p.Pmin, p.k);
    dSSF_Px4 = aux_fcn_dSSF(xk4(4), p.Pmax, p.Pmin, p.k);
    Idpre    = SSF_Px4 / Vf;

    dIdpre_dx1 = -SSF_Px4 / Vf^2;
    dIdpre_dx4 = dSSF_Px4 / Vf;

    SSF_Id     = aux_fcn_SSF(Idpre, p.Idmax, p.Idmin, p.k);
    dSSF_Id    = aux_fcn_dSSF(Idpre, p.Idmax, p.Idmin, p.k);

    Jf4(5,1) = dt * dSSF_Id * dIdpre_dx1 / Tg;
    Jf4(5,4) = dt * dSSF_Id * dIdpre_dx4 / Tg;
    Jf4(5,5) = -dt / Tg;

    % % === Optional Logging ===
    % fprintf('[Jf4] k=%d | SDBF=%.3e | dIqv_dkqv=%.3e | SSF_P=%.3e | SSF_Id=%.3e\n', ...
    %         k, SDBF, dIqv_dkqv, SSF_P, SSF_Id);
end
