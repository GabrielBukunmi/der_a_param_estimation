function k4 = compute_k4_vradynamics(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now)
    xk4 = xk + k3;
    k4 = compute_k1_vradynamics(xk4, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
end
