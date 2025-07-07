function [Vmult] = pcode_trip_volt(p,x)
Vt = x(1);
Vmin = Vt;
Timer1 = 0;
Timer2 = 0;
Counter1 = 0;
Counter2 = 0;
if Vt < p.vl1 & Timer1 == 0
    Timer1 = tic;
elseif Vt > p.vl1 & Timer1 > 0
    Timer1 = 0;
end
if Vt < p.vl0 & Timer2 == 0
    Timer2 = tic;
elseif Vt > p.vl0 & Timer2 > 0
    Timer2 = 0;
end
if Vmin <= p.vl0
    Vmin = p.vl0;
end
if Vt <= p.vl0 | Counter2 == 1
    Multiplier = 0;
elseif Vt <= p.vl1 & Counter1 == 0
    Multiplier = (Vt-p.vl0) / (p.vl1-p.vl0);
elseif Vt <= p.vl1 & Counter1 == 1
    Multiplier = ((Vmin-p.vl0) +p.Vrfrac*(Vt-Vmin)) / (p.vl1-p.vl0);
elseif Vt >= p.vl1 & Counter1 == 0
    Multiplier = 1;
else
    Multiplier = p.Vrfrac*((p.vl1-Vmin)/(p.vl1-p.vl0)) + ((Vmin-p.vl0)/(p.vl1-p.vl0));
end
if Counter1 == 0
    if Timer1 > p.tvl1
        Counter1 = 1;
        p.Vmin = Vt;
    end
end
if Counter2 == 0
    if Timer2 > p.tvl0
        Counter2 = 1;
    end
end
Vmult = Vt * Multiplier;