function k4 = compute_k4_P_f_control(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    % === Full-step update ===
    xk4 = xk + k3;

    % === Evaluate k1 at full-step ===
    k4 = compute_k1_P_f_control(xk4, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
end
