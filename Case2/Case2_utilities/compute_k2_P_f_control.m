function k2 = compute_k2_P_f_control(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    % === Midpoint update ===
    xk2 = xk + 0.5 * k1;

    % === Evaluate k1 at midpoint ===
    k2 = compute_k1_P_f_control(xk2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
end
