function [A] = aux_fcn_A_2(p, x, Pref_now)
    A = aux_fcn_SSF((Pref_now + aux_fcn_uD(p, x) - x(2)), p.femax, p.femin, p.k);
end
