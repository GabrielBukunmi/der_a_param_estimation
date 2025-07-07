function [B] = aux_fcn_B(p,x)
B = aux_fcn_SSF( x(7) + p.kpg*aux_fcn_A(p,x),p.Pmax,p.Pmin,p.k);
end