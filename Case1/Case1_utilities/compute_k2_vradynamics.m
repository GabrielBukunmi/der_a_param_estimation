function k2 = compute_k2_vradynamics(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    xk2 = xk + 0.5 * k1;  % full update: both states and parameters
    k2 = compute_k1_vradynamics(xk2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
end
