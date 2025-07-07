function Idpre = aux_fcn_idpre_p_f_control(p, x_Idpre)
    % Compute Idpre based on:
    % x_Idpre(1): Vt_filt (filtered terminal voltage)
    % x_Idpre(2): Pgen (active power generation)

    eps_denom = 1e-4;

    Vt_filt = max(x_Idpre(1), eps_denom);      % Avoid divide-by-zero
    Pgen    = x_Idpre(2);                      % Active power state

    Pclip   = aux_fcn_SSF(Pgen, p.Pmax, p.Pmin, p.k);  % Saturated Pgen
    Idpre   = Pclip / Vt_filt;                         % Compute P/V
end
