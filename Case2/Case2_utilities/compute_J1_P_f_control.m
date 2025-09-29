function J1 = compute_J1_P_f_control(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    % Jacobian matrix for stage 1 of RK4 integration (J1)
    % 14 states: [x1..x8] dynamics, [x9..x14] parameters (held constant)

    J1 = zeros(14, 14);
    eps = 1e-6;
    idx = max(k-1, 1);

    % Inputs
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    freq  = freq_now(idx);

    % States
    Vtfilt    = max(xk(1), eps);
    PgenFilt  = xk(2);
    IqPFC     = xk(3);
    Iq        = xk(4);
    FreqFilt  = xk(5);
    PI_Int    = xk(6);
    Pgen      = xk(7);
    Id        = xk(8);

    % Parameters (from state vector)
    Tp   = max(xk(9), eps);
    Kpg  = xk(10);
    Kig  = xk(11);
    Trf  = max(xk(12), eps);
    Ddn  = xk(13);
    Dup  = xk(14);

    %% Row 1: d(Vtfilt)/dt
    J1(1,1) = -dt / p.Trv;

    %% Row 2: d(PgenFilt)/dt
    dSSF_Pgen = aux_fcn_dSSF(Pgen, p.Pmax, p.Pmin, p.k);
    J1(2,2) = -dt / Tp;
    J1(2,7) = dt * dSSF_Pgen / Tp;
    J1(2,9) = dt * (PgenFilt - aux_fcn_SSF(Pgen, p.Pmax, p.Pmin, p.k)) / Tp^2;

    %% Row 3: d(IqPFC)/dt
    pfaref = atan(Qref / (Pref + eps));
    Vclip = aux_fcn_SSF(Vtfilt, p.inf, 0.01, p.k);
    dVclip = aux_fcn_dSSF(Vtfilt, p.inf, 0.01, p.k);
    J1(3,1) = -dt * pfaref * PgenFilt * dVclip / (p.Tiq * Vclip^2);
    J1(3,2) = dt * pfaref / (p.Tiq * Vclip);
    J1(3,3) = -dt / p.Tiq;

    %% Row 4: d(Iq)/dt
    Iqv    = aux_fcn_Iqv(p, xk);
    Iqref  = aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k);
    dIqref = aux_fcn_dSSF(Iqv, p.Iqh1, p.Iql1, p.k);
    psi    = IqPFC - Iqref;
    dIqcmd = aux_fcn_dSSF(psi, p.Iqmax, p.Iqmin, p.k);
    J1(4,3) = dt * dIqcmd / p.Tg;
    J1(4,4) = -dt / p.Tg;

    %% Row 5: d(FreqFilt)/dt
    J1(5,5)  = -dt / Trf;
    J1(5,12) = dt * (FreqFilt - freq) / Trf^2;

    %% Row 6: d(PI_Integral)/dt
    Ferr    = p.Freq_ref - FreqFilt;
    SDBF1   = aux_fcn_SDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
    dSDBF1  = -aux_fcn_dSDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
    dSSF_dn = aux_fcn_dSSF(SDBF1 * Ddn, 0, -p.inf, p.k);
    dSSF_up = aux_fcn_dSSF(SDBF1 * Dup, p.inf, 0, p.k);

    uD = aux_fcn_SSF(SDBF1 * Ddn, 0, -p.inf, p.k) + aux_fcn_SSF(SDBF1 * Dup, p.inf, 0, p.k);
    A = aux_fcn_SSF(Pref + uD - PgenFilt, p.Pmax, p.Pmin, p.k);
    B = aux_fcn_SSF(PI_Int + Kpg * A, p.Pmax, p.Pmin, p.k);

    dA = aux_fcn_dSSF(Pref + uD - PgenFilt, p.Pmax, p.Pmin, p.k);
    dB = aux_fcn_dSSF(PI_Int + Kpg * A, p.Pmax, p.Pmin, p.k);

    J1(6,2)  = -dt * Kig * dA - dt * p.kw * dB * (-Kpg * dA);
    J1(6,6)  = -dt * p.kw * dB;
    J1(6,10) = -dt * p.kw * dB * A;
    J1(6,11) = dt * dA;
    J1(6,5)  = dt * (-Ddn * dSSF_dn * dSDBF1 - Dup * dSSF_up * dSDBF1);
    J1(6,13) = dt * dSSF_dn * SDBF1;
    J1(6,14) = dt * dSSF_up * SDBF1;

    %% Row 7: d(Pgen)/dt
    Pclip   = aux_fcn_SSF(Pgen, p.Pmax, p.Pmin, p.k);
    dPclip  = aux_fcn_dSSF(Pgen, p.Pmax, p.Pmin, p.k);
    J1(7,7)  = dt * (-dPclip / p.Tpord - p.kw);

    %% Row 8: d(Id)/dt
    Vsafe    = max(Vtfilt, eps);
    Pclip_id = aux_fcn_SSF(Pgen, p.Pmax, p.Pmin, p.k);
    Idpre    = Pclip_id / Vsafe;
    dIdpre_V = -Pclip_id / (Vsafe^2);
    dIdpre_P = aux_fcn_dSSF(Pgen, p.Pmax, p.Pmin, p.k) / Vsafe;
    dIdcmd   = aux_fcn_dSSF(Idpre, p.Idmax, p.Idmin, p.k);

    J1(8,1) = dt * dIdcmd * dIdpre_V / p.Tg;
    J1(8,7) = dt * dIdcmd * dIdpre_P / p.Tg;
    J1(8,8) = -dt / p.Tg;

    % Rows 9–14 (parameters) have ∂dot/∂param = 0 (frozen), so remain zero
end
