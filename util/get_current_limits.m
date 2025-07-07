function [Ipmax, Ipmin, Iqmax, Iqmin] = get_current_limits(p, Ipcmd, Iqcmd)
    if p.PQflag == 0  % Q-priority
        Iqmax = p.Imax;
        Iqmin = -p.Imax;
        Ipmax = sqrt(max(p.Imax^2 - Iqcmd^2, 0));
    else  % P-priority
        Ipmax = p.Imax;
        Iqmax = sqrt(max(p.Imax^2 - Ipcmd^2, 0));
        Iqmin = -Iqmax;
    end
    Ipmin = (p.typeflag == 0) * (-Ipmax);
end
