function J1 = compute_J1_vradynamics(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    J1 = zeros(10, 10);
    idx = max(k-1, 1);
    eps_denom = 1e-6;

    % === Safe parameter extraction ===
    Trv   = max(xk(6), eps_denom);
    kqv   = max(xk(7), eps_denom);
    Tg    = max(xk(8), eps_denom);
    Tiq   = max(xk(9), eps_denom);
    Tpord = max(xk(10), eps_denom);

    Vt   = Vt_now(idx);
    Qref = Qref_now(idx);
    Pref = Pref_now(idx);
    Vf   = max(xk(1), eps_denom);

    % === Row 1: d(x1)/dt ===
    J1(1,1) = -dt / Trv;
    J1(1,6) = dt * (xk(1) - Vt) / Trv^2;

    % === Row 2: d(iqpfc)/dt ===
    J1(2,1) = dt * Qref / (Tiq * Vf^2);
    J1(2,2) = -dt / Tiq;
    J1(2,9) = dt * (xk(2) - Qref / Vf) / Tiq^2;

    % === Row 3: d(x3)/dt === (reactive control)
    deltaV     = p.Vref - xk(1);
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
    psi        = xk(2) - SSF_Iqv;
    dSSF_psi   = aux_fcn_dSSF(psi, p.Iqmax, p.Iqmin, p.k);

    J1(3,1) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dV) / Tg;
    J1(3,2) = dt * dSSF_psi / Tg;
    J1(3,3) = -dt / Tg;
    J1(3,7) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dkqv) / Tg;
    J1(3,8) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dbd2) / Tg;
    % Optional: J1(3,9) if dbd1 is estimated too:
    % J1(3,9) = dt * dSSF_psi * (-dSSF_Iqv * dIqv_dbd1) / Tg;

    % === Row 4: d(x4)/dt ===
    J1(4,4)  = -dt / Tpord;
    J1(4,10) = dt * (xk(4) - Pref) / Tpord^2;

    % === Row 5: d(x5)/dt ===
    J1(5,5) = -dt / Tg;

    % % === Debug Print ===
    % fprintf('[J1] k=%d | dIqv_dkqv=%.2e | dSDBF_dbd2=%.2e | dSSF_Iqv=%.2e\n', ...
    %         k, dIqv_dkqv, dSDBF_dbd2, dSSF_Iqv);
end
