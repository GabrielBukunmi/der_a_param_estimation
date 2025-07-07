function k1 = compute_k1_P_f_control(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)

    % === Initialize ===
    k1 = zeros(14,1);  % 8 states + 6 parameters
    eps = 1e-6;
    idx = max(k-1, 1);

    % === Clamp ===
    xk = min(max(xk, -1e8), 1e6);
    xk(9:14) = max(xk(9:14), eps);  % Only 6 parameters

    % === Inputs ===
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    Pref  = Pref_now(idx);
    freq  = freq_now(idx);

    % === Extract Parameters ===
    Tp    = xk(9);
    Kpg   = xk(10);
    Kig   = xk(11);
    Trf   = xk(12);
    Ddn   = xk(13);
    Dup   = xk(14);

    % === x1: VtFilt ===
    k1(1) = dt * (-1/p.Trv * xk(1) + 1/p.Trv * Vt);

    % === x2: PgenFilt ===
    try
        k1(2) = dt * (-1/Tp * xk(2) + 1/Tp * aux_fcn_SSF(xk(7), p.Pmax, p.Pmin, p.k));
    catch
        k1(2) = 0;
        warning('[k1] PgenFilt SSF failed at k=%d.', k);
    end

    % === x3: IqPFC ===
    try
        pfaref = atan(Qref / (Pref + eps));
        Vt_clip = aux_fcn_SSF(xk(1), p.inf, 0.01, p.k);
        k1(3) = dt * (-1/p.Tiq * xk(3) + 1/p.Tiq * pfaref * xk(2) / Vt_clip);
    catch
        k1(3) = 0;
        warning('[k1] IqPFC failed at k=%d.', k);
    end

    % === x4: Iq ===
    try
        Iqv    = aux_fcn_Iqv(p, xk);
        Iqref  = aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k);
        Iqcmd  = aux_fcn_SSF(xk(3) - Iqref, p.Iqmax, p.Iqmin, p.k);
        k1(4) = dt * (-1/p.Tg * xk(4) + 1/p.Tg * Iqcmd);
    catch
        k1(4) = 0;
        warning('[k1] Iq control failed at k=%d.', k);
    end

    % === x5: FreqFilt ===
    k1(5) = dt * (-1/Trf * xk(5) + 1/Trf * freq);

    % === x6: PI_Integral ===
    try
        Ferr = p.Freq_ref - xk(5);
        SDBF1 = aux_fcn_SDBF(Ferr, p.fdbd2, p.fdbd1, p.k);
        uD = aux_fcn_SSF(SDBF1 * Ddn, 0, -p.inf, p.k) + aux_fcn_SSF(SDBF1 * Dup, p.inf, 0, p.k);

        A = aux_fcn_SSF(Pref + uD - xk(2), p.femax, p.femin, p.k);
        B = aux_fcn_SSF(xk(6) + Kpg * A, p.Pmax, p.Pmin, p.k);
        k1(6) = dt * (Kig * A + p.kw * (B - Kpg * A - xk(6)));
    catch
        k1(6) = 0;
        warning('[k1] PI Integral (A/B/uD) failed at k=%d.', k);
    end

    % === x7: Pgen ===
    try
        Pclip = aux_fcn_SSF(xk(7), p.Pmax, p.Pmin, p.k);
        k1(7) = dt * (1/p.Tpord * B - 1/p.Tpord * Pclip + p.kw * (Pclip - xk(7)));
    catch
        k1(7) = 0;
        warning('[k1] Pgen update failed at k=%d.', k);
    end

    % === x8: Id ===
    try
        x_Idpre = [xk(1); xk(7)];  % x1: Vt_filt, x7: Pgen
        Idpre   = aux_fcn_idpre_p_f_control(p, x_Idpre);
        Idcmd   = aux_fcn_SSF(Idpre, p.Idmax, p.Idmin, p.k);
        k1(8)   = dt * (1/p.Tg * Idcmd - 1/p.Tg * xk(8));
    catch
        k1(8) = 0;
        warning('[k1] Id control failed at k=%d.', k);
    end

    % === Parameters: Frozen ===
    k1(9:14) = 0;

    % === Final Clamp ===
    k1 = min(max(k1, -1e6), 1e6);
end
