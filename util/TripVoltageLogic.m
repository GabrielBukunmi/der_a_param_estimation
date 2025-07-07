function Vmult = TripVoltageLogic(Time, Vt, V10, V11, Vh0, Vh1, Vrfrac, Tvl1, Tvl0, Tvh1, Tvh0)
    % Persistent state
    persistent Vmin Vmax TimeBelowV11 TimeBelowV10 TimeAboveVh1 TimeAboveVh0
    persistent ActiveTimerV11 ActiveTimerV10 ActiveTimerVh1 ActiveTimerVh0
    persistent ActiveFracLow ActiveFracHigh ActiveTripLow ActiveTripHigh

    % Initialize on first run
    if isempty(Vmin)
        Vmin = V11;
        Vmax = Vh1;
        TimeBelowV11 = 0; TimeBelowV10 = 0;
        TimeAboveVh1 = 0; TimeAboveVh0 = 0;
        ActiveTimerV11 = false; ActiveTimerV10 = false;
        ActiveTimerVh1 = false; ActiveTimerVh0 = false;
        ActiveFracLow = false; ActiveFracHigh = false;
        ActiveTripLow = false; ActiveTripHigh = false;
    end

    %% ----- LOW VOLTAGE RIDE-THROUGH -----
    if Vt >= V11
        ActiveTimerV11 = false;
    elseif ~ActiveTimerV11
        ActiveTimerV11 = true;
        TimeBelowV11 = Time;
    end

    if Vt >= V10
        ActiveTimerV10 = false;
    elseif ~ActiveTimerV10
        ActiveTimerV10 = true;
        TimeBelowV10 = Time;
    end

    if ~ActiveFracLow && ActiveTimerV11 && (Time - TimeBelowV11) >= Tvl1
        ActiveFracLow = true;
    end

    if ~ActiveTripLow && ActiveTimerV10 && (Time - TimeBelowV10) >= Tvl0
        ActiveTripLow = true;
    end

    if Vmin > Vt && ActiveFracLow
        Vmin = Vt;
    end

    if Vt <= V10 || ActiveTripLow
        VmultL = 0;
    elseif Vt <= V11 && ActiveFracLow && Vt > Vmin
        VmultL = ((Vt - V10) + Vrfrac * (Vmin - V10)) / (V11 - V10);
    elseif Vt < V11
        VmultL = (Vt - V10) / (V11 - V10);
    elseif ~ActiveFracLow
        VmultL = 1.0;
    else
        VmultL = ((Vt - V10) + Vrfrac * (V11 - Vmin)) / (V11 - V10);
    end

    %% ----- HIGH VOLTAGE RIDE-THROUGH -----
    if Vt <= Vh1
        ActiveTimerVh1 = false;
    elseif ~ActiveTimerVh1
        ActiveTimerVh1 = true;
        TimeAboveVh1 = Time;
    end

    if Vt <= Vh0
        ActiveTimerVh0 = false;
    elseif ~ActiveTimerVh0
        ActiveTimerVh0 = true;
        TimeAboveVh0 = Time;
    end

    if ~ActiveFracHigh && ActiveTimerVh1 && (Time - TimeAboveVh1) >= Tvh1
        ActiveFracHigh = true;
    end

    if ~ActiveTripHigh && ActiveTimerVh0 && (Time - TimeAboveVh0) >= Tvh0
        ActiveTripHigh = true;
    end

    if Vmax < Vt && ActiveFracHigh
        Vmax = Vt;
    end

    if Vt >= Vh0 || ActiveTripHigh
        VmultH = 0;
    elseif Vt >= Vh1 && ActiveFracHigh && Vt < Vmax
        VmultH = ((Vt - Vh0) + Vrfrac * (Vh1 - Vmax)) / (Vh1 - Vh0);
    elseif Vt < Vh1
        VmultH = (Vt - Vh0) / (Vh1 - Vh0);
    elseif ~ActiveFracHigh
        VmultH = 1.0;
    else
        VmultH = ((Vt - Vh0) + Vrfrac * (Vh1 - Vmax)) / (Vh1 - Vh0);
    end

    %% ----- FINAL MULTIPLIER -----
    Vmult = min(VmultL, VmultH);  % DER follows stricter limit (lower multiplier)
end
