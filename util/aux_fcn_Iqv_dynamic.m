function [Iqv] = aux_fcn_Iqv_dynamic(Vref_now, Vt_now, dbd2, dbd1, kqv, k)
    deadband_out = aux_fcn_SDBF(Vref_now - Vt_now, dbd2, dbd1, k);
    Iqv = kqv * deadband_out;
end
