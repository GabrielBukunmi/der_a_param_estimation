function Jf3 = compute_Jf3_P_f_control(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    Jf3 = zeros(14, 14);
    idx = max(k-1, 1);
    eps = 1e-6;

    % RK4 midpoint state update
    xk3 = xk;
    xk3(1:8) = xk(1:8) + 0.5 * k2(1:8);

    % Inputs
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    freq  = freq_now(idx);

    % Parameters (ensure positive)
    Tp   = max(xk3(9), eps);
    Kpg  = max(xk3(10), eps);
    Kig  = max(xk3(11), eps);
    Trf  = max(xk3(12), eps);
    Ddn  = xk3(13);
    Dup  = xk3(14);

    %% Row 1: VtFilt
    Jf3(1,1) = -dt / p.Trv;

    %% Row 2: PgenFilt
    dSSF_Pgen = aux_fcn_dSSF(xk3(7), p.Pmax, p.Pmin, p.k);
    Jf3(2,2) = -dt / Tp;
    Jf3(2,7) = dt * dSSF_Pgen / Tp;
    Jf3(2,9) = dt * (xk3(2) - aux_fcn_SSF(xk3(7), p.Pmax, p.Pmin, p.k)) / Tp^2;

    %% Row 3: IqPFC
    pfaref = atan(Qref / (Pref + eps));
    Vt_clip = aux_fcn_SSF(xk3(1), p.inf, 0.01, p.k);
    dSSF_Vt = aux_fcn_dSSF(xk3(1), p.inf, 0.01, p.k);
    Jf3(3,1) = -dt * pfaref * xk3(2) * dSSF_Vt / (p.Tiq * Vt_clip^2);
    Jf3(3,2) = dt * pfaref / (p.Tiq * Vt_clip);
    Jf3(3,3) = -dt / p.Tiq;

    %% Row 4: Iq
    Iqv     = aux_fcn_Iqv(p, xk3);
    Iqref   = aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k);
    dIqref  = aux_fcn_dSSF(Iqv, p.Iqh1, p.Iql1, p.k);
    psi     = xk3(3) - Iqref;
    dpsi    = aux_fcn_dSSF(psi, p.Iqmax, p.Iqmin, p.k);
    Jf3(4,3) = dt * dpsi / p.Tg;
    Jf3(4,4) = -dt / p.Tg;

    %% Row 5: FreqFilt
    Jf3(5,5) = -dt / Trf;
    Jf3(5,12) = dt * (xk3(5) - freq) / Trf^2;

    %% Row 6: PI_Integral
    A = aux_fcn_A_2(p, xk3, Pref);
    B = aux_fcn_B(p, xk3);
    dA_d2 = -aux_fcn_dSSF(xk3(2), p.Pmax, p.Pmin, p.k);
    dB_d7 = aux_fcn_dSSF(xk3(7) + Kpg * A, p.Pmax, p.Pmin, p.k);
    dB_A  = aux_fcn_dSSF(B, p.Pmax, p.Pmin, p.k);

    Jf3(6,2) = dt * (Kig * dA_d2 - p.kw * Kpg * dA_d2);
    Jf3(6,6) = -dt * p.kw;
    Jf3(6,7) = dt * p.kw * dB_d7;
    Jf3(6,10) = dt * p.kw * A * dB_A;
    Jf3(6,11) = dt * A;

    % uD derivatives
    Ferr     = p.Freq_ref - xk3(5);
    SDBF1    = aux_fcn_SDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
    dSDBF1   = -aux_fcn_dSDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
    dSSF_dn  = aux_fcn_dSSF(SDBF1 * Ddn, 0, -p.inf, p.k);
    dSSF_up  = aux_fcn_dSSF(SDBF1 * Dup, p.inf, 0, p.k);

    Jf3(6,5)  = dt * (-Ddn * dSSF_dn * dSDBF1 - Dup * dSSF_up * dSDBF1);
    Jf3(6,13) = dt * dSSF_dn * SDBF1;
    Jf3(6,14) = dt * dSSF_up * SDBF1;

    %% Row 7: Pgen
    Pclip = aux_fcn_SSF(xk3(7), p.Pmax, p.Pmin, p.k);
    dPclip = aux_fcn_dSSF(xk3(7), p.Pmax, p.Pmin, p.k);
    Jf3(7,7) = dt * (-dPclip / p.Tpord - p.kw);

    %% Row 8: Id
    Vtfilt_safe = max(xk3(1), eps);
    Pclip_id = aux_fcn_SSF(xk3(7), p.Pmax, p.Pmin, p.k);
    Idpre = Pclip_id / Vtfilt_safe;
    dIdpre_dVt = -Pclip_id / (Vtfilt_safe^2);
    dIdpre_dP  = aux_fcn_dSSF(xk3(7), p.Pmax, p.Pmin, p.k) / Vtfilt_safe;
    dSSF_Id = aux_fcn_dSSF(Idpre, p.Idmax, p.Idmin, p.k);

    Jf3(8,1) = dt * dSSF_Id * dIdpre_dVt / p.Tg;
    Jf3(8,7) = dt * dSSF_Id * dIdpre_dP / p.Tg;
    Jf3(8,8) = -dt / p.Tg;

    %% Rows 9–14: (parameters) remain zero
    % ∂param_dot/∂xk = 0 for frozen parameters
end
