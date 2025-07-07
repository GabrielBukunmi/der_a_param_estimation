function Jf2 = compute_Jf2_P_f_control(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    % Jacobian for midpoint RK4 (k2) for P-F control: 8 states + 6 parameters

    Jf2 = zeros(14, 14);
    eps = 1e-6;
    idx = max(k-1, 1);

    % Midpoint state estimate
    xk2 = xk;
    xk2(1:8) = xk(1:8) + 0.5 * k1(1:8);

    % Extract states
    Vtfilt    = max(xk2(1), eps);
    PgenFilt  = xk2(2);
    IqPFC     = xk2(3);
    Iq        = xk2(4);
    FreqFilt  = xk2(5);
    PI_Int    = xk2(6);
    Pgen      = xk2(7);
    Id        = xk2(8);

    % Extract parameters
    Tp   = max(xk2(9), eps);
    Kpg  = xk2(10);
    Kig  = xk2(11);
    Trf  = max(xk2(12), eps);
    Ddn  = xk2(13);
    Dup  = xk2(14);

    % Inputs
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    freq  = freq_now(idx);

    %% Row 1: VtFilt
    Jf2(1,1) = -dt / p.Trv;

    %% Row 2: PgenFilt
    dSSF_Pgen = aux_fcn_dSSF(Pgen, p.Pmax, p.Pmin, p.k);
    Jf2(2,2) = -dt / Tp;
    Jf2(2,7) = dt * dSSF_Pgen / Tp;
    Jf2(2,9) = dt * (PgenFilt - aux_fcn_SSF(Pgen, p.Pmax, p.Pmin, p.k)) / Tp^2;

    %% Row 3: IqPFC
    pfaref   = atan(Qref / (Pref + eps));
    V_clip   = aux_fcn_SSF(Vtfilt, p.inf, 0.01, p.k);
    dSSF_V   = aux_fcn_dSSF(Vtfilt, p.inf, 0.01, p.k);
    Jf2(3,1) = -dt * pfaref * PgenFilt * dSSF_V / (p.Tiq * V_clip^2);
    Jf2(3,2) = dt * pfaref / (p.Tiq * V_clip);
    Jf2(3,3) = -dt / p.Tiq;

    %% Row 4: Iq
    Iqv        = aux_fcn_Iqv(p, xk2);
    dIqref     = aux_fcn_dSSF(Iqv, p.Iqh1, p.Iql1, p.k);
    dIqcmd     = aux_fcn_dSSF(IqPFC - aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k), p.Iqmax, p.Iqmin, p.k);
    Jf2(4,3)   = dt * dIqcmd / p.Tg;
    Jf2(4,4)   = -dt / p.Tg;

    %% Row 5: FreqFilt
    Jf2(5,5)    = -dt / Trf;
    Jf2(5,12)   = dt * (FreqFilt - freq) / Trf^2;

    %% Row 6: PI_Integral
    A      = aux_fcn_A_2(p, xk2, Pref);
    B      = aux_fcn_B(p, xk2);
    dA_d2  = -aux_fcn_dSSF(xk2(2), p.Pmax, p.Pmin, p.k);
    dB_d7  = aux_fcn_dSSF(xk2(7) + Kpg * A, p.Pmax, p.Pmin, p.k);
    dB     = aux_fcn_dSSF(B, p.Pmax, p.Pmin, p.k);

    Jf2(6,2)  = dt * (Kig * dA_d2 - p.kw * Kpg * dA_d2);
    Jf2(6,6)  = -dt * p.kw;
    Jf2(6,7)  = dt * p.kw * dB_d7;
    Jf2(6,10) = dt * p.kw * A * dB;
    Jf2(6,11) = dt * A;

    % Partial wrt frequency-related droop control
    Ferr      = p.Freq_ref - FreqFilt;
    dSDBF     = aux_fcn_dSDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
    sign_dn   = double(Ferr < 0);
    sign_up   = double(Ferr > 0);
    SDBF_val  = aux_fcn_SDBF(Ferr, p.fdbd2, p.fdbd1, p.k);

    Jf2(6,5)   = dt * (-Ddn * dSDBF * sign_dn - Dup * dSDBF * sign_up);
    Jf2(6,13)  = dt * SDBF_val * sign_dn;
    Jf2(6,14)  = dt * SDBF_val * sign_up;

    %% Row 7: Pgen
    dSSF_Pgen = aux_fcn_dSSF(Pgen, p.Pmax, p.Pmin, p.k);
    Jf2(7,7)  = dt * (-1 / p.Tpord - p.kw) * dSSF_Pgen;

    %% Row 8: Id
    Idpre      = aux_fcn_Idpre_vra(p, [Vtfilt; Pgen]);
    dIdpre_dV  = -aux_fcn_SSF(Pgen, p.Pmax, p.Pmin, p.k) / (Vtfilt^2);
    dIdpre_dP  = aux_fcn_dSSF(Pgen, p.Pmax, p.Pmin, p.k) / Vtfilt;
    dIdcmd     = aux_fcn_dSSF(Idpre, p.Idmax, p.Idmin, p.k);
    Jf2(8,1)   = dt * dIdcmd * dIdpre_dV / p.Tg;
    Jf2(8,7)   = dt * dIdcmd * dIdpre_dP / p.Tg;
    Jf2(8,8)   = -dt / p.Tg;

    %% Rows 9–14 (parameters): zero rows
    % Parameters are frozen → ∂x_dot/∂param = 0
end
