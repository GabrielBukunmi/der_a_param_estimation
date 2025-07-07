function Idpre = aux_fcn_Idpre_vra(p, x)
    % x(1): voltage (x1), x(2): power control state (x4)
    eps_denom = 1e-4;
    V = max(x(1),eps_denom);
    P_ctrl = aux_fcn_SSF(x(2), p.Pmax, p.Pmin, p.k);

    Idpre = P_ctrl / V;
end

