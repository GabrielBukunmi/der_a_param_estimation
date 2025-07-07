function [uD] = aux_fcn_uD(p, x)
    Ferr = p.Freq_ref - x(6);  % Frequency error

    term1 = aux_fcn_SDBF(Ferr, p.fdbd2, p.fdbd1, p.k) * p.Ddn;
    term2 = aux_fcn_SDBF(Ferr, p.fdbd2, p.fdbd1, p.k) * p.Dup;

    uD = aux_fcn_SSF(term1, 0, -p.inf, p.k) + aux_fcn_SSF(term2, p.inf, 0, p.k);
end
