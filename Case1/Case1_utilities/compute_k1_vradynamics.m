function k1 = compute_k1_vradynamics(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)

    % === Initialize ===
    k1 = zeros(10,1);  % 5 states + 5 parameters
    eps_denom = 1e-6;
    idx = max(k-1, 1);

    % === Clamp full xk to avoid NaN/Inf ===
    xk = min(max(xk, -1e6), 1e6);    % hard bound on all states/params
    xk(6:10) = max(xk(6:10), eps_denom);  % ensure all parameters positive

    % === Extract Parameters from Augmented State ===
    Trv   = xk(6);
    kqv   = xk(7);
    Tg    = xk(8);
    Tiq   = xk(9);
    Tpord = xk(10);

    % === Inputs at Previous Time Step ===
    Vt    = Vt_now(idx);
    Qref  = Qref_now(idx);
    freq  = freq_now(idx);  % Unused but reserved
    Pref  = Pref_now(idx);

    % === Voltage Filter (x1) ===
    Vfilter = max(xk(1), eps_denom);  % avoid divide-by-zero

    % === Reactive Current Control ===
    x_v = xk;
    p.kqv = kqv;  % override p.kqv with current estimate

    % -- Guard for Iqv and derivatives inside helper functions
    try
        Iqv     = aux_fcn_Iqv(p, x_v);  % kqv * SDBF
        Iq_ref  = aux_fcn_SSF(Iqv, p.Iqh1, p.Iql1, p.k);
        Iq_sat  = aux_fcn_SSF(xk(2) - Iq_ref, p.Iqmax, p.Iqmin, p.k);
    catch
        warning('[k1] SDF or SSF failed at k=%d. Resetting Iqv.', k);
        Iqv = 0; Iq_ref = 0; Iq_sat = 0;
    end

    % === D-axis Current Control ===
    x_Idpre = [xk(1); xk(4)];  % x1: V_filtered, x4: Pref control state
    try
        Idpre   = aux_fcn_Idpre_vra(p, x_Idpre);
        Id_sat  = aux_fcn_SSF(Idpre, p.Idmax, p.Idmin, p.k);
    catch
        warning('[k1] Idpre or SSF failed at k=%d. Resetting Id_sat.', k);
        Id_sat = 0;
    end

    % === Dynamic State Derivatives ===
    k1(1) = dt * (-xk(1)/Trv + Vt/Trv);                         % x1: V_filtered
    k1(2) = dt * (-xk(2)/Tiq + Qref / (Tiq * Vfilter));         % x2: iqpfc
    k1(3) = dt * (-xk(3)/Tg + Iq_sat/Tg);                       % x3: Iq

    % Pref controller (x4)
    try
        Pctrl     = aux_fcn_SSF(xk(4), p.Pmax, p.Pmin, p.k);
        dPref     = Pref / Tpord - Pctrl / Tpord + p.kw * (Pctrl - xk(4));
    catch
        warning('[k1] Pref controller SSF failed at k=%d.', k);
        dPref = 0;
    end
    k1(4) = dt * dPref;                                        % x4: Pref

    % Id control (x5)
    k1(5) = dt * (Id_sat / Tg - xk(5)/Tg);                     % x5: Id

    % === Parameter Derivatives (set to 0) ===
    k1(6:10) = 0;

    % === Clamp Output ===
    k1 = min(max(k1, -1e6), 1e6);  % avoid extreme steps
end
