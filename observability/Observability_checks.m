clear; clc;

addpath(genpath('../util'))

load('der_a_DATA.mat');

t = Vmag_800.Time;
[t_unique, unique_idx] = unique(t, 'stable');
t = t_unique;

Vmag_802_pu = Vmag_802.Data(unique_idx) / ((24900/sqrt(3)) * sqrt(2));
Q_total = (Q_PV_824.Data(unique_idx) + Q_PV_842.Data(unique_idx) + ...
           Q_PV_850.Data(unique_idx) + Q_PV_832.Data(unique_idx) + ...
           Q_PV_840.Data(unique_idx)) / 100;
P_total = (P_PV_824.Data(unique_idx) + P_PV_842.Data(unique_idx) + ...
           P_PV_850.Data(unique_idx) + P_PV_832.Data(unique_idx) + ...
           P_PV_840.Data(unique_idx)) / 100;

Vt_now_raw = Vmag_802_pu(1);
Vt_now = max(Vt_now_raw, 1e-3);
Qref_now = Q_total(1);
Pref_now = P_total(1);

%% 
nx = 6;
ntheta = 17;
n = nx + ntheta;

syms x [n 1] real

% States
x1 = x(1); x3 = x(2); x4 = x(3); x8 = x(4); x9 = x(5); x10 = x(6);
% Parameters
Trv = x(7); Tpord = x(8); Tiq = x(9); Tg = x(10);
kqv = x(11); Iqh1 = x(12); Iql1 = x(13);
Pmax = x(14); Pmin = x(15); Idmax = x(16); Idmin = x(17);
dpmax = x(18); dpmin = x(19); dbd1 = x(20); dbd2 = x(21);
Iqmax = x(22); Iqmin = x(23);

% Constants
Ts = 0.033;
kw = 10;

%% 
% Smooth Saturation Function
SSF = @(z, U, L, k) ...
    ((U + L)/2) + ((U - L)/2) * ( (z - (U + L)/2) / ((U - L)/2 + 1e-6) ) ./ ...
    ((1 + abs((z - (U + L)/2) / ((U - L)/2 + 1e-6)).^k).^(1/k));

% Smooth Deadband Function
SDBF = @(z, U, L, k) ...
    (z - (U + L)/2) - ((U - L)/2) * ( (z - (U + L)/2) / ((U - L)/2 + 1e-6) ) ./ ...
    ((1 + abs((z - (U + L)/2) / ((U - L)/2 + 1e-6)).^k).^(1/k));

k= 1024;

%% === Model Equations ===

Iqv = kqv * SDBF(x1 - 0.9, dbd2, dbd1, k);
Iq_inner = SSF(Iqv, Iqh1, Iql1, k);
Iq_eff = SSF(x3 - Iq_inner, Iqmax, Iqmin, k);
Pref_raw = -x8 / Ts + Pref_now / Ts;
Pref_filt = SSF(Pref_raw, dpmax, dpmin, k);
Pref_sat = SSF(x9, Pmax, Pmin, k);
Id_pre = Pref_sat / x1;
Id_eff = SSF(Id_pre, Idmax, Idmin, k);

%%
f = sym(zeros(n,1));
f(1) = -x1 / Trv + Vt_now / Trv;
f(2) = -x3 / Tiq + Qref_now / (Tiq * x1);
f(3) = -x4 / Tg + Iq_eff / Tg;
f(4) = Pref_filt;
f(5) = (x8 - Pref_sat)/Tpord + kw * (Pref_sat - x9);
f(6) = (Id_eff - x10) / Tg;
f(7:end) = 0;
%% 
outputs = {x1; x1*x10; x1*x4};

%% 
O = [];
for j = 1:length(outputs)
    h = outputs{j};
    Lf0 = simplify(jacobian(h, x));
    O   = [O; Lf0];
    for k = 1:1
        Lf_val = simplify(Lf0 * f);
        Lf0 = simplify(jacobian(Lf_val, x));
        O = [O; Lf0];
    end
end

%% Display and Save
rank_O = rank(O);
disp(['Rank of Observability Matrix: ', num2str(rank_O)]);

save('observability_symbolic_only.mat', 'O', 'rank_O');
