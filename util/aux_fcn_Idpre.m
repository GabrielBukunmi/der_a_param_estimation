function [Idpre] = aux_fcn_Idpre(p,x)
Idpre = aux_fcn_SSF(x(9), p.Pmax,p.Pmin,p.k) / x(1);
end