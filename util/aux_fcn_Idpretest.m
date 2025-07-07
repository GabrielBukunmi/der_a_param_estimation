function [Idpre] = aux_fcn_Idpretest(p,x)
Idpre = aux_fcn_SSF(x(5), p.Pmax,p.Pmin,p.k) / aux_fcn_SSF(x(1), p.inf,0.01,p.k);
end