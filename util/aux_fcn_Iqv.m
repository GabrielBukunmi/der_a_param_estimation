function [Iqv] = aux_fcn_Iqv(p, x)
    deadband_out = aux_fcn_SDBF(p.Vref - x(1), p.dbd2, p.dbd1, p.k);
    Iqv = p.kqv * deadband_out;
end
