function k3 = compute_k3_vradynamics(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    xk3 = xk + 0.5 * k2;
    k3 = compute_k1_vradynamics(xk3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
end
