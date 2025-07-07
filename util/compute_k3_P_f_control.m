function k3 = compute_k3_P_f_control(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    % === Midpoint update ===
    xk3 = xk + 0.5 * k2;

    % === Evaluate k1 at midpoint ===
    k3 = compute_k1_P_f_control(xk3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
end
