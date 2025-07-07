function Jf4 = compute_Jf4_P_f_control(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    Jf4 = zeros(14, 14);
    idx = max(k-1, 1);
    eps = 1e-6;

    % RK4 final stage state update
    xk4 = xk;
    xk4(1:8) = xk(1:8) + k3(1:8);

    % Inputs
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    freq  = freq_now(idx);

    % Parameters (ensure positive where needed)
    Tp   = max(xk4(9), eps);
    Kpg  = max(xk4(10), eps);
    Kig  = max(xk4(11), eps);
    Trf  = max(xk4(12), eps);
    Ddn  = xk4(13);
    Dup  = xk4(14);

    %% Row 1: VtFilt
    Jf4(1,1) = -dt / p.Trv;

    %% Row 2: PgenFilt
    dSSF_Pgen = aux_fcn_dSSF(xk4(7), p.Pmax, p.Pmin, p.k);
    Jf4(2,2) = -dt / Tp;
    Jf4(2,7) = dt * dSSF_Pgen / Tp;
    Jf4(2,9) = dt * (xk4(2) - aux_fcn_SSF(xk4(7), p.Pmax, p.Pmin, p.k)) / Tp^2;

    %% Row 3: IqPFC
    pfaref = atan(Qref / (Pref + eps));
    Vt_clip = aux_fcn_SSF(xk4(1), p.inf, 0.01, p.k);
    dSSF_Vt = aux_fcn_dSSF(xk4(1), p.inf, 0.01, p.k);
    Jf4(3,1) = -dt * pfaref * xk4(2) * dSSF_Vt / (p.Tiq * Vt_clip^2);
    Jf4(3,2) = dt * pfaref / (p.Tiq * Vt_clip);
    Jf4(3,3) = -dt / p.Tiq;

    %% Row 4: Iq
    Iqv     = aux_fcn_Iqv(p, xk4);
    Iqref   = aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k);
    dIqref  = aux_fcn_dSSF(Iqv, p.Iqh1, p.Iql1, p.k);
    psi     = xk4(3) - Iqref;
    dpsi    = aux_fcn_dSSF(psi, p.Iqmax, p.Iqmin, p.k);
    Jf4(4,3) = dt * dpsi / p.Tg;
    Jf4(4,4) = -dt / p.Tg;

    %% Row 5: FreqFilt
    Jf4(5,5) = -dt / Trf;
    Jf4(5,12) = dt * (xk4(5) - freq) / Trf^2;

    %% Row 6: PI_Integral
    A = aux_fcn_A_2(p, xk4, Pref);
    B = aux_fcn_B(p, xk4);
    dA_d2 = -aux_fcn_dSSF(xk4(2), p.Pmax, p.Pmin, p.k);
    dB_d7 = aux_fcn_dSSF(xk4(7) + Kpg * A, p.Pmax, p.Pmin, p.k);
    dB_A  = aux_fcn_dSSF(B, p.Pmax, p.Pmin, p.k);

    Jf4(6,2) = dt * (Kig * dA_d2 - p.kw * Kpg * dA_d2);
    Jf4(6,6) = -dt * p.kw;
    Jf4(6,7) = dt * p.kw * dB_d7;
    Jf4(6,10) = dt * p.kw * A * dB_A;
    Jf4(6,11) = dt * A;

    Ferr     = p.Freq_ref - xk4(5);
    SDBF1    = aux_fcn_SDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
    dSDBF1   = -aux_fcn_dSDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
    dSSF_dn  = aux_fcn_dSSF(SDBF1 * Ddn, 0, -p.inf, p.k);
    dSSF_up  = aux_fcn_dSSF(SDBF1 * Dup, p.inf, 0, p.k);

    Jf4(6,5)  = dt * (-Ddn * dSSF_dn * dSDBF1 - Dup * dSSF_up * dSDBF1);
    Jf4(6,13) = dt * dSSF_dn * SDBF1;
    Jf4(6,14) = dt * dSSF_up * SDBF1;

    %% Row 7: Pgen
    Pclip = aux_fcn_SSF(xk4(7), p.Pmax, p.Pmin, p.k);
    dPclip = aux_fcn_dSSF(xk4(7), p.Pmax, p.Pmin, p.k);
    Jf4(7,7) = dt * (-dPclip / p.Tpord - p.kw);

    %% Row 8: Id
    Vtfilt_safe = max(xk4(1), eps);
    Pclip_id = aux_fcn_SSF(xk4(7), p.Pmax, p.Pmin, p.k);
    Idpre = Pclip_id / Vtfilt_safe;
    dIdpre_dVt = -Pclip_id / (Vtfilt_safe^2);
    dIdpre_dP  = aux_fcn_dSSF(xk4(7), p.Pmax, p.Pmin, p.k) / Vtfilt_safe;
    dSSF_Id = aux_fcn_dSSF(Idpre, p.Idmax, p.Idmin, p.k);

    Jf4(8,1) = dt * dSSF_Id * dIdpre_dVt / p.Tg;
    Jf4(8,7) = dt * dSSF_Id * dIdpre_dP / p.Tg;
    Jf4(8,8) = -dt / p.Tg;

    %% Done
end
