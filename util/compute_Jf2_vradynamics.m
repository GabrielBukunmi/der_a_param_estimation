function Jf2 = compute_Jf2_vradynamics(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    Jf2 = zeros(10, 10);
    eps_denom = 1e-6;
    idx = max(k-1, 1);

    % RK4 midpoint state update
    xk2 = xk;
    xk2(1:5) = xk(1:5) + 0.5 * k1(1:5);

    % Parameter extraction
    Trv   = max(xk2(6), eps_denom);
    kqv   = max(xk2(7), eps_denom);
    Tg    = max(xk2(8), eps_denom);
    Tiq   = max(xk2(9), eps_denom);
    Tpord = max(xk2(10), eps_denom);

    % Inputs
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    Vf    = max(xk2(1), eps_denom);

    % === Row 1: d(x1)/dt ===
    Jf2(1,1) = -dt / Trv;
    Jf2(1,6) = dt * (xk2(1) - Vt) / Trv^2;

    % === Row 2: d(iqpfc)/dt ===
    Jf2(2,1) = dt * Qref / (Tiq * Vf^2);
    Jf2(2,2) = -dt / Tiq;
    Jf2(2,9) = dt * (xk2(2) - Qref / Vf) / Tiq^2;

    % === Row 3: d(Iq)/dt (reactive dynamics) ===
    deltaV     = p.Vref - xk2(1);
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
    psi        = xk2(2) - SSF_Iqv;
    dSSF_psi   = aux_fcn_dSSF(psi, p.Iqmax, p.Iqmin, p.k);

    Jf2(3,1) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dV) / Tg;
    Jf2(3,2) = dt * dSSF_psi / Tg;
    Jf2(3,3) = -dt / Tg;
    Jf2(3,7) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dkqv) / Tg;
    Jf2(3,8) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dbd2) / Tg;
    % Optional: Jf2(3,9) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dbd1) / Tg;

    % === Row 4: d(Pord)/dt ===
    Jf2(4,4)  = -dt / Tpord;
    Jf2(4,10) = dt * (xk2(4) - Pref) / Tpord^2;

    % === Row 5: d(Id)/dt ===
    Jf2(5,5) = -dt / Tg;

    % % === Debugging ===
    % fprintf('[Jf2] k = %d | dIqv_dkqv = %.2e | dSSF_Iqv = %.2e\n', ...
    %         k, dIqv_dkqv, dSSF_Iqv);
end
