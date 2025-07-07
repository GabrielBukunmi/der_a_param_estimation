function [A] = aux_fcn_A(p,x)
A = aux_fcn_SSF((p.Pref + aux_fcn_uD(p,x)-x(2)), p.femax,p.femin,p.k);
end