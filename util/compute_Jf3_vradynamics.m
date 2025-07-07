function Jf3 = compute_Jf3_vradynamics(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    Jf3 = zeros(10, 10);
    idx = max(k-1, 1);
    eps_denom = 1e-6;

    % Half-step state update
    xk3 = xk;
    xk3(1:5) = xk(1:5) + 0.5 * k2(1:5);

    % Parameters
    Trv   = max(xk(6), eps_denom);
    kqv   = max(xk(7), eps_denom);
    Tg    = max(xk(8), eps_denom);
    Tiq   = max(xk(9), eps_denom);
    Tpord = max(xk(10), eps_denom);

    % Inputs
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    Vf    = max(xk3(1), eps_denom);

    % === Row 1: d(s0)/dt ===
    Jf3(1,1) = -dt / Trv;
    Jf3(1,6) = dt * (xk3(1) - Vt) / Trv^2;

    % === Row 2: d(iqpfc)/dt ===
    Jf3(2,1) = dt * Qref / (Tiq * Vf^2);
    Jf3(2,2) = -dt / Tiq;
    Jf3(2,9) = dt * (xk3(2) - Qref / Vf) / Tiq^2;

    % === Row 3: d(Iq)/dt (reactive control) ===
    deltaV     = p.Vref - xk3(1);
    dbd1       = max(p.dbd1, eps_denom);
    dbd2       = max(p.dbd2, eps_denom);

    SDBF       = aux_fcn_SDBF(deltaV, dbd2, dbd1, p.k);
    dSDBF_dV   = -aux_fcn_dSDBF(deltaV, dbd2, dbd1, p.k);
    dSDBF_dbd2 = aux_fcn_dSDBF_dbd2(deltaV, dbd2, dbd1, p.k);
    dSDBF_dbd1 = aux_fcn_dSDBF_dbd1(deltaV, dbd2, dbd1, p.k);

    Iqv        = kqv * SDBF;
    dIqv_dV    = kqv * dSDBF_dV;
    dIqv_dkqv  = SDBF;
    dIqv_dbd2  = kqv * dSDBF_dbd2;
    dIqv_dbd1  = kqv * dSDBF_dbd1;

    SSF_Iqv    = aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k);
    dSSF_Iqv   = aux_fcn_dSSF(Iqv, p.Iqh1, p.Iql1, p.k);
    psi        = xk3(2) - SSF_Iqv;
    dSSF_psi   = aux_fcn_dSSF(psi, p.Iqmax, p.Iqmin, p.k);

    Jf3(3,1) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dV) / Tg;
    Jf3(3,2) = dt * dSSF_psi / Tg;
    Jf3(3,3) = -dt / Tg;
    Jf3(3,7) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dkqv) / Tg;
    Jf3(3,8) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dbd2) / Tg;
    % Optional: Jf3(3,9) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dbd1) / Tg;

    % === Row 4: d(Pord)/dt ===
    Jf3(4,4)  = -dt / Tpord;
    Jf3(4,10) = dt * (xk3(4) - Pref) / Tpord^2;

    % === Row 5: d(Id)/dt ===
    Jf3(5,5) = -dt / Tg;

    % % === Logging ===
    % fprintf('[Jf3] k=%d | Iqv=%.2e | dIqv_dV=%.2e | dIqv_dbd2=%.2e\n', ...
    %     k, Iqv, dIqv_dV, dIqv_dbd2);
end
