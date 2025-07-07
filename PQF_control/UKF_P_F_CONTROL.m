%% Pflag =1 Vtrfl=0
% Freqflag =1

clear; close all; clc
addpath(genpath('../util'))

%% 

 load('der_a_DATA.mat')
t = Vmag_800.Time;

% === Step 2: Remove duplicate timestamps ===
[t_unique, unique_idx] = unique(t, 'stable');

% === Step 3: Apply unique indices to all signals ===
Vmag_800 = Vmag_800.Data(unique_idx);
Vmag_802 = Vmag_802.Data(unique_idx);
Imag_800 = Imag_800.Data(unique_idx);
Imag_802 = Imag_802.Data(unique_idx);
Vang_800 = Vang_800.Data(unique_idx);
Vang_802 = Vang_802.Data(unique_idx);
Iang_800 = Iang_800.Data(unique_idx);
Iang_802 = Iang_802.Data(unique_idx);

Vdq0_800 = Vdq0_800.Data(unique_idx, :);
Vdq0_802 = Vdq0_802.Data(unique_idx, :);
Idq0_800 = Idq0_800.Data(unique_idx, :);
Idq0_802 = Idq0_802.Data(unique_idx, :);
f_802 = f_802.Data(unique_idx, :);
P_800 = P_800.Data(unique_idx);
Q_800 = Q_800.Data(unique_idx);
Pv_P = Pv_P.Data(unique_idx);
Qv_P = Qv_P.Data(unique_idx);

P_PV_842 = P_PV_842.Data(unique_idx);
P_PV_850 = P_PV_850.Data(unique_idx);
P_PV_832 = P_PV_832.Data(unique_idx);
P_PV_840 = P_PV_840.Data(unique_idx);
P_PV_824 = P_PV_824.Data(unique_idx);

Q_PV_824 = Q_PV_824.Data(unique_idx);
Q_PV_840 = Q_PV_840.Data(unique_idx);
Q_PV_832 = Q_PV_832.Data(unique_idx);
Q_PV_842 = Q_PV_842.Data(unique_idx);
Q_PV_850 = Q_PV_850.Data(unique_idx);

t = t_unique;  % Assign cleaned time vector

% === Step 4: Now trim from t >= 0.2 seconds ===
start_idx = find(t >= 0.05, 1);

t = t(start_idx:end);
Vmag_800 = Vmag_800(start_idx:end);
Vmag_802 = Vmag_802(start_idx:end);
Imag_800 = Imag_800(start_idx:end);
Imag_802 = Imag_802(start_idx:end);
Vang_800 = Vang_800(start_idx:end);
Vang_802 = Vang_802(start_idx:end);
Iang_800 = Iang_800(start_idx:end);
Iang_802 = Iang_802(start_idx:end);
Vdq0_800 = Vdq0_800(start_idx:end, :);
Vdq0_802 = Vdq0_802(start_idx:end, :);
Idq0_800 = Idq0_800(start_idx:end, :);
Idq0_802 = Idq0_802(start_idx:end, :);
f_802 = f_802(start_idx:end, :);
P_800 = P_800(start_idx:end);
Q_800 = Q_800(start_idx:end);
Pv_P = Pv_P(start_idx:end);
Qv_P = Qv_P(start_idx:end);
%% PV DATA

P_PV_842 = P_PV_842(start_idx:end);
P_PV_850 = P_PV_850(start_idx:end);
P_PV_832 = P_PV_832(start_idx:end);
P_PV_840 = P_PV_840(start_idx:end);
P_PV_824 = P_PV_824(start_idx:end);

Q_PV_824 = Q_PV_824(start_idx:end);
Q_PV_840 = Q_PV_840(start_idx:end);
Q_PV_832 = Q_PV_832(start_idx:end);
Q_PV_842 = Q_PV_842(start_idx:end);
Q_PV_850 = Q_PV_850(start_idx:end);

%% 
% per unit base values

Fbase =60;
Vbase = 24900/ sqrt(3) * sqrt(2);   % Peak phase voltage from 25 kV L-L RMS
Sbase = 100e3;                      % 100 kVA
Ibase= Sbase/ Vbase;
%% 

% positive sequence data in per unit and radians
Vmag_800_pu = Vmag_800 ./ Vbase;
Vmag_802_pu = Vmag_802 ./ Vbase;
Imag_800_pu = Imag_800 ./ Ibase;
Imag_802_pu = Imag_802 ./ Ibase;

Vang_800_rad = Vang_800 * pi/180;
Vang_802_rad = Vang_802 * pi/180;
Iang_800_rad = Iang_800 * pi/180;
Iang_802_rad = Iang_802 * pi/180;

% % positive sequence plots
% figure % plot 1
% subplot(2,2,1); plot(t, Vmag_800); hold on; plot(t, Vmag_802); legend('Vmag 800 (kV)','Vmag 802 (kV)')
% subplot(2,2,2); plot(t, Vang_800); hold on; plot(t, Vang_802); legend('Vang 800 (deg)','Vang 802 (deg)')
% subplot(2,2,3); plot(t, Vmag_800 - Vmag_802); legend('Subtraction: Vmag 800 - Vmag 802 (kV)')
% subplot(2,2,4); plot(t, Vang_800 - Vang_802); legend('Subtraction: Vang 800 - Vang 802 (deg)')
% 
% figure % plot 2
% subplot(2,2,1); plot(t, Vmag_800_pu);  hold on; plot(t, Vmag_802_pu);  legend('Vmag 800 (pu)','Vmag 802 (pu)')
% subplot(2,2,2); plot(t, Vang_800_rad); hold on; plot(t, Vang_802_rad); legend('Vang 800 (rad)','Vang 802 (rad)')
% subplot(2,2,3); plot(t, Vmag_800_pu  - Vmag_802_pu);  legend('Subtraction: Vmag 800 - Vmag 802 (pu)')
% subplot(2,2,4); plot(t, Vang_800_rad - Vang_802_rad); legend('Subtraction: Vang 800 - Vang 802 (rad)')
% 
% figure % plot 3
% subplot(2,2,1); plot(t, Imag_800); hold on; plot(t, Imag_802); legend('Imag 800 (A)','Imag 802 (A)')
% subplot(2,2,2); plot(t, Iang_800); hold on; plot(t, Iang_802); legend('Iang 800 (deg)','Iang 802 (deg)')
% subplot(2,2,3); plot(t, Imag_800 - Imag_802); legend('Subtraction: Imag 800 - Imag 802 (A)')
% subplot(2,2,4); plot(t, Iang_800 - Iang_802); legend('Subtraction: Iang 800 - Iang 802 (deg)')
% 
% figure % plot 4
% subplot(2,2,1); plot(t, Imag_800_pu);  hold on; plot(t, Imag_802_pu);  legend('Imag 800 (pu)','Imag 802 (pu)')
% subplot(2,2,2); plot(t, Iang_800_rad); hold on; plot(t, Iang_802_rad); legend('Iang 800 (rad)','Iang 802 (rad)')
% subplot(2,2,3); plot(t, Imag_800_pu  - Imag_802_pu);  legend('Subtraction: Imag 800 - Imag 802 (pu)')
% subplot(2,2,4); plot(t, Iang_800_rad - Iang_802_rad); legend('Subtraction: Iang 800 - Iang 802 (rad)')

% dq0 data in per unit
Vd_800_pu = Vdq0_800(:,1) ./ Vbase;
Vq_800_pu = Vdq0_800(:,2) ./ Vbase;
V0_800_pu = Vdq0_800(:,3) ./ Vbase;

Vd_802_pu = Vdq0_802(:,1) ./ Vbase;
Vq_802_pu = Vdq0_802(:,2) ./ Vbase;
V0_802_pu = Vdq0_802(:,3) ./ Vbase;

Id_800_pu = Idq0_800(:,1) ./ Ibase;
Iq_800_pu = Idq0_800(:,2) ./ Ibase;
I0_800_pu = Idq0_800(:,3) ./ Ibase;

Id_802_pu = Idq0_802(:,1) ./ Ibase;
Iq_802_pu = Idq0_802(:,2) ./ Ibase;
I0_802_pu = Idq0_802(:,3) ./ Ibase;
% 
% % dq0 plots
% figure % plot 5
% subplot(2,3,1); plot(t, Vd_800_pu); legend('Vd 800 (pu)')
% subplot(2,3,2); plot(t, Vq_800_pu); legend('Vq 800 (pu)')
% subplot(2,3,3); plot(t, V0_800_pu); legend('V0 800 (pu)')
% subplot(2,3,4); plot(t, Vd_802_pu); legend('Vd 802 (pu)')
% subplot(2,3,5); plot(t, Vq_802_pu); legend('Vq 802 (pu)')
% subplot(2,3,6); plot(t, V0_802_pu); legend('V0 802 (pu)')

% dq0 voltages and currents
V_800_pu  = Vd_800_pu + 1i * Vq_800_pu;
V_800_mag = abs(V_800_pu);
V_800_ang = angle(V_800_pu);

V_802_pu  = Vd_802_pu + 1i * Vq_802_pu;
V_802_mag = abs(V_802_pu);
V_802_ang = angle(V_802_pu);

I_800_pu  = Id_800_pu + 1i * Iq_800_pu;
I_800_mag = abs(I_800_pu);
I_800_ang = angle(I_800_pu);

I_802_pu  = Id_802_pu + 1i * Iq_802_pu;
I_802_mag = abs(I_802_pu);
I_802_ang = angle(I_802_pu);

% dq0 vs. positive sequence plots
% figure % plot 6
% subplot(2,3,1); plot(t, V_800_mag); hold on; plot(t, Vmag_800_pu); legend('Vmag800 dq0 (pu)','Vmag800 pos. seq. (pu)')
% subplot(2,3,2); plot(t, V_800_mag - Vmag_800_pu); legend('Subtraction: Vmag800 dq0 - Vmag800 pos. seq. (pu)')
% subplot(2,3,3); plot(t, V_802_mag - Vmag_802_pu); legend('Subtraction: Vmag802 dq0 - Vmag802 pos. seq. (pu)')
% subplot(2,3,4); plot(t, I_800_mag); hold on; plot(t, Imag_800_pu); legend('Imag800 dq0 (pu)','Imag800 pos. seq. (pu)')
% subplot(2,3,5); plot(t, I_800_mag - Imag_800_pu); legend('Subtraction: Imag800 dq0 - Imag800 pos. seq. (pu)')
% subplot(2,3,6); plot(t, I_802_mag - Imag_802_pu); legend('Subtraction: Imag802 dq0 - Imag802 pos. seq. (pu)')
% 
% figure % plot 7
% subplot(2,3,1); plot(t, V_800_ang); hold on; plot(t, Vang_800_rad); legend('Vang800 dq0 (rad)','Vang800 pos. seq. (rad)')
% subplot(2,3,2); plot(t, V_800_ang - Vang_800_rad); legend('Subtraction: Vang800 dq0 - Vang800 pos. seq. (rad)')
% subplot(2,3,3); plot(t, V_802_ang - Vang_802_rad); legend('Subtraction: Vang802 dq0 - Vang802 pos. seq. (rad)')
% subplot(2,3,4); plot(t, I_800_ang); hold on; plot(t, Iang_800_rad); legend('Iang800 dq0 (rad)','Iang800 pos. seq. (rad)')
% subplot(2,3,5); plot(t, I_800_ang - Iang_800_rad); legend('Subtraction: Iang800 dq0 - Iang800 pos. seq. (rad)')
% subplot(2,3,6); plot(t, I_802_ang - Iang_802_rad); legend('Subtraction: Iang802 dq0 - Iang802 pos. seq. (rad)')

% Power
P800pu_simu = P_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
Q800pu_simu = Q_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
P800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* cos(Vang_800_rad - Iang_800_rad);
Q800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* sin(Vang_800_rad - Iang_800_rad);

% Compute active and reactive power (pu) at node 802
P802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* cos(Vang_802_rad - Iang_802_rad);
Q802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* sin(Vang_802_rad - Iang_802_rad);


% % Power plots
% figure % plot 8
% subplot(1,2,1); plot(t, P800pu_simu,'k-'); hold on; plot(t, P800pu_dq0,'r--'); legend('P from Simulink pos. seq. block (pu)','P calculated in MATLAB pos. seq. (pu)')
% subplot(1,2,2); plot(t, Q800pu_simu,'k-'); hold on; plot(t, Q800pu_dq0,'r--'); legend('Q from Simulink pos. seq. block (pu)','Q calculated in MATLAB pos. seq. (pu)')

%% 

% Sum
P_total = (P_PV_824 + P_PV_842 + P_PV_850 + P_PV_832 + P_PV_840)/100; %Already dividing the data by 1000 in simulink
Q_total = (Q_PV_824 + Q_PV_842 + Q_PV_850 + Q_PV_832 + Q_PV_840)/100;





% 
% figure % plot 3 – total PV active and reactive power
% subplot(2,1,1); plot(t, P_total, 'LineWidth', 1.2); title('Total PV Active Power'); ylabel('P_{total} (kW)'); grid on
% subplot(2,1,2); plot(t, Q_total, 'LineWidth', 1.2); title('Total PV Reactive Power'); ylabel('Q_{total} (kVAR)'); grid on
% xlabel('Time (s)');



% Define DER-A parameters 
p = struct(...                   % Reference: 2023 - NERC - Parameterization of the DER_A Model for Aggregate DER.pdf
    'dbd1'     , -0.05     , ... % lower voltage deadband ≤ 0 (pu)
    'dbd2'     ,  0.05     , ... % upper voltage deadband ≥ 0 (pu)
    'Ddn'      , 20      , ... % frequency control droop gain ≥ 0 (down-side)
    'dpmax'    ,  25      , ... % power ramp rate up > 0 (pu/s)
    'dpmin'    , -25      , ... % power ramp rate down < 0 (pu/s)
    'Dup'      ,  20      , ... % frequency control droop gain ≥ 0 (up-side)
    'fdbd1'    , -0.000002 , ... % lower frequency control deadband ≤ 0 (pu)
    'fdbd2'    ,  0.000002 , ... % upper frequency control deadband ≥ 0 (pu)
    'femax'    ,  99       , ... % frequency control maximum error ≥ 0 (pu)
    'femin'    , -99       , ... % frequency control minimum error ≤ 0 (pu)
    'Fh'       ,  60.07    , ... % frequency break-point for high frequency cut-out of inverters
    'Fl'       ,  59.93    , ... % frequency break-point for low frequency cut-out of inverters
    'Freq'     ,  1        , ... % frequency
    'Freq_ref' ,  1        , ... % set-point frequency reference
    'inf'      ,  99    , ... % infinity (algebraic variable)
    'neg_inf'  , -99   ,.... % negative infinity (algebraic variable)
    'Idmax'    ,  10.9      , ... % maximum limit of total active current (or Ipmax)
    'Idmin'    , -10.9      , ... % minimum limit of total active current (or Ipmin)
    'Imax'     ,  10.2      , ... % maximum converter current (pu)
    'Iqh1'     ,  2.2        , ... % maximum limit of reactive current injection, p.u.
    'Iql1'     , -2.2        , ... % minimum limit of reactive current injection, p.u.    
    'Iqmax'    ,  4.0      , ... % maximum limit of total reactive current
    'Iqmin'    , -4.0     , ... % minimum limit of total reactive current
    'k'        ,  1024     , ... % smoothing factor
    'kig'      ,  10       , ... % active power control integral gain
    'kpg'      ,  0.5      , ... % active power control proportional gain
    'kqv'      ,  5       , ... % proportional voltage control gain (pu/pu)
    'kw'       ,  2      , ... % time constant for anti-windup (algbraic)
    'Pmax'     ,   8       , ... % maximum power (pu)
    'Pmin'     , 0        , ... % minimum power (pu)
    'Pref'     ,  0.9      , ... % active power reference at point of fault
    'Qref'     ,  0.1   , ... % reactive power reference
    'rrpwr'    ,  0.05    , ... % Power rise ramp rate following a fault > 0 (pu/s)
    'rrdn'     ,  7        , ... % (created by the paper) - Power rise ramp rate down side, if s9>=0, rrdn = -∞, otherwise, rrdn= -rrpwr
    'rrup'     ,  7       , ... % (created by the paper) - Power rise ramp rate up side, if s9>=0, rrdn = rrpw, otherwise, rrup= ∞
    'Tfh'      ,  4.5      , ... % timer for fh
    'Tfl'      ,  4.6      , ... % timer for fl (Tfl > Trf)
    'Tg'       ,  0.02     , ... % current control time constant
    'Tiq'      ,  0.02     , ... % Q control time constant (s)
    'Tp'       ,  0.02    , ... % transducer time constant (s)
    'Tpord'    ,  0.02     , ... % power order time constant (s)
    'Trf'      ,  0.02     , ... % transducer time constant (seconds) for frequency measurement (> 0)  
    'Trv'      ,  0.02     , ... % transducer time constant (s) for voltage measurement
    'Ts'       ,  0.02     , ... % Evaluation time of input signal (created by the author of the paper for simplification purposes)
    'Tv'       ,  0.02     , ... % time constant on the output of the voltage/frequency cut-out
    'tvh0'     ,  1.5     , ... % timer for vh0 point
    'tvh1'     ,  1.5        , ... % timer for vh1 point
    'tvl0'     ,  1.5     , ... % timer for vl0 point
    'tvl1'     ,  2        , ... % timer for vl1 point
    'vh0'      ,  1      , ... % voltage break-point for high voltage cut-out of inverters
    'vh1'      ,  1.05     , ... % voltage break-point for high voltage cut-out of inverters
    'vl0'      ,  0.8      , ... % voltage break-point for low voltage cut-out of inverters
    'vl1'      ,  0.91     , ... % voltage break-point for low voltage cut-out of inverters
    'Vpr'      ,  0.9      , ... % voltage below which frequency tripping is disabled
    'Vref'     ,  0.9     , ... % voltage reference set-point > 0 (pu)
    'Vrfrac'   ,  0.6      , ... % fraction of device that recovers after voltage comes back to within vl1 < V < vh1
    'Xe'       ,  0.009      ,  ... % source impedance reactive > 0 (pu)
    'pfaref'   ,  0.87   ,    ... % power factor reference for when Pflag=1
    'PQflag'    , 0 ,....    % 0 for Q-priority and 1 for P priority
    'typeflag' , 0 .....      %TypeFlag: 0 means the unit is a storage device and Ipmin = - Ipmax; 1 means the unit is a generator Ipmin = 0 (any number which is not 0 is treated as 1)
    );


% Voltage and current at t = t(1)
Vpu   = Vmag_802_pu(1);
Ipu   = Imag_802_pu(1);
Vang  = Vang_802_rad(1);
Iang  = Iang_802_rad(1);


%% % State initialization


% State initialization
s0_0 = Vpu;
s7_0 = P_total(1);
s8_0 = solveForS8_0(p, s7_0);
s1_0 = aux_fcn_SSF(s8_0, p.Pmax, p.Pmin, p.k);
s2_0 = Q_total(1);
Iqv_0 = aux_fcn_Iqv(p, s0_0);                      % From deadband-based function
Iq_ref_0 = aux_fcn_SSF(Iqv_0, p.Iqh1, p.Iql1, p.k); % Saturated reference
s3_0 = aux_fcn_SSF(s2_0 - Iq_ref_0, p.Iqmax, p.Iqmin, p.k);
s5_0 = f_802(1) / Fbase;

% The initialization of s6 requires an iterative numerical method
x(2) = s1_0; 
x(6) = s5_0; 
x(7) = 0.2; % Initial guess for s6_0
x_updated = solveForS6_0(p, x);
s6_0 = x_updated(7);


x(1) = s0_0;  
x(9) = s8_0;
s9_0 = 0.5;


%% True state propagation with augmented parameter 
%% Power-Frequency Control: Initialization (8 States + 6 Parameters)

n_states = 8;         % Dynamic states: x1 to x8
n_param  = 6;         % Parameters to estimate: x9 to x14
n_total  = n_states + n_param;

dt = 0.02;            % Time step (s)
Nt = length(t);       % Number of simulation steps

% === Preallocate State Matrix ===
x_true = zeros(n_total, Nt);

% === Initial Conditions (States + Parameters) ===
x_true(:,1) = [ ...
    s0_0;         % x1: Vt_filt
    s1_0;         % x2: PgenFilt
    s2_0;         % x3: IqPFC
    s3_0;         % x4: Iq
    s5_0;         % x5: FreqFilt
    s6_0;         % x6: PI_Integral
    s7_0;         % x7: Pgen
    s9_0;         % x8: Id
    p.Tp;         % x9: Active power filter time constant
    p.kpg;        % x10: Frequency loop proportional gain
    p.kig;        % x11: Frequency loop integral gain
    p.Trf;        % x12: Frequency filter time constant
    p.Ddn;        % x13: Downward droop gain
    p.Dup         % x14: Upward droop gain
];

% === Inputs ===
Vt_now   = Vmag_802_pu;       % Terminal voltage (per unit)
Pref_now = P_total;           % Active power reference
Qref_now = Q_total;           % Reactive power reference
freq_now = f_802 / 60;        % Normalized frequency

%% RK4 Propagation of Dynamics
for k = 2:Nt-1
    xk = x_true(:, k-1);

    k1 = compute_k1_P_f_control(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k2 = compute_k2_P_f_control(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k3 = compute_k3_P_f_control(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k4 = compute_k4_P_f_control(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);

    x_true(:, k) = x_true(:, k-1) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end

%% Plotting Trajectories
tt = t(1:Nt-1);  % Trimmed time vector for plotting

state_labels = {
    'x1: Vt_filt';
    'x2: PgenFilt';
    'x3: IqPFC';
    'x4: Iq';
    'x5: FreqFilt';
    'x6: PI_Integral';
    'x7: Pgen';
    'x8: Id';
    'x9: Tp';
    'x10: Kpg';
    'x11: Kig';
    'x12: Trf';
    'x13: Ddn';
    'x14: Dup'
};

figure;
for i = 1:n_total
    subplot(4, 4, i);
    plot(tt, x_true(i, 1:Nt-1), 'LineWidth', 1.5);
    grid on;
    title(state_labels{i}, 'Interpreter', 'none');
    xlabel('Time (s)');
    ylabel(state_labels{i});
end

sgtitle('RK4 Trajectories of Power-Frequency States and Parameters');


%% %% === UKF Initialization ===
n = 14;   % Total: 8 states + 6 parameters
alpha = 1e-3; kappa = 0; beta = 2;
lambda = alpha^2 * (n + kappa) - n;
gamma = sqrt(n + lambda);

% UKF weights
Wm = [lambda / (n + lambda), repmat(1 / (2 * (n + lambda)), 1, 2*n)];
Wc = Wm; Wc(1) = Wc(1) + (1 - alpha^2 + beta);

% Preallocation
xhat_UKF = zeros(n, Nt);
P_UKF = zeros(n, n, Nt);

% Initial state
xhat_UKF(:,1) = [s0_0; s1_0; s2_0; s3_0; s5_0; s6_0; s7_0; s9_0; ...
                 p.Tp - 0.01; p.kpg * 0.95; p.kig - 0.03; ...
                 p.Trf - 0.0005; p.Ddn - 0.5; p.Dup * 0.97];

% Initial covariance
P_UKF(:,:,1) = diag([0.01 * ones(1,8), 0.1 * ones(1,6)]);

% Process and measurement noise
Q = (1e-3)^2 * eye(n);
R = (1e-4)^2 * eye(6);  % 6 measurements

%% === UKF Estimation Loop ===
for k = 2:Nt-1
    xk = xhat_UKF(:,k-1);
    Pk = P_UKF(:,:,k-1);

    % --- Sigma Points ---
    A = gamma * chol(Pk + 1e-9*eye(n), 'lower');
    X_sigma = [xk, xk + A, xk - A];   % 2n+1 columns

    % --- Propagate each sigma point using RK4 ---
    X_sigma_pred = zeros(n, 2*n+1);
    for i = 1:(2*n+1)
        xi = X_sigma(:,i);

        k1 = compute_k1_P_f_control(xi, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
        k2 = compute_k2_P_f_control(xi, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
        k3 = compute_k3_P_f_control(xi, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
        k4 = compute_k4_P_f_control(xi, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);

        X_sigma_pred(:,i) = xi + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    % --- Predicted Mean and Covariance ---
    x_pred = X_sigma_pred * Wm';
    P_pred = Q;
    for i = 1:(2*n+1)
        dx = X_sigma_pred(:,i) - x_pred;
        P_pred = P_pred + Wc(i) * (dx * dx');
    end

    % --- Predict Measurement for Each Sigma Point ---
    Y_sigma = zeros(6, 2*n+1);
    for i = 1:(2*n+1)
        xi = X_sigma_pred(:,i);
        x1 = xi(1); x4 = xi(4); x5 = xi(5); x6 = xi(6); x7 = xi(7);
        Tp   = xi(9);
        Kpg  = xi(10);
        Kig  = xi(11);
        Trf  = xi(12);
        Ddn  = xi(13);
        Dup  = xi(14);

        Ferr = p.Freq_ref - x5;
        uD = Ferr * Ddn + Ferr * Dup;
        dy5 = (-1/Trf) * x5 + (1/Trf) * freq_now(k);

        Y_sigma(:,i) = [
            x1;                     % Vt
            x1 * x4;                % Q = Vt * Iq
            x5;                     % Frequency
            x7;                     % Pgen
            dy5;                    % df/dt
            uD                     
        ];
    end
    y_pred = Y_sigma * Wm';

    % --- Covariances ---
    Pyy = R;
    Pxy = zeros(n, 6);
    for i = 1:(2*n+1)
        dy = Y_sigma(:,i) - y_pred;
        dx = X_sigma_pred(:,i) - x_pred;
        Pyy = Pyy + Wc(i) * (dy * dy');
        Pxy = Pxy + Wc(i) * (dx * dy');
    end

    % --- Kalman Gain and Measurement ---
    Kk = Pxy / Pyy;
    x_true_k = x_true(:,k);
    y_meas = [
        x_true_k(1);
        x_true_k(1) * x_true_k(4);
        x_true_k(5);
        x_true_k(7);
        (x_true_k(5) - x_true(5,k-1)) / dt;
        (p.Freq_ref - x_true_k(5)) * x_true_k(13) + ...
        (p.Freq_ref - x_true_k(5)) * x_true_k(14)
    ];

    % --- Update ---
    xhat_UKF(:,k) = x_pred + Kk * (y_meas - y_pred);
    xhat_UKF(1:8,k) = min(max(xhat_UKF(1:8,k), -1e6), 1e6);     % Clamp states
    xhat_UKF(9:end,k) = max(xhat_UKF(9:end,k), 1e-6);           % Parameters positive

    P_UKF(:,:,k) = P_pred - Kk * Pyy * Kk';
end
% === UKF Estimation-Error Variance Tracking ===
ukf_var = zeros(6, Nt);
for i = 1:6
    for k = 1:Nt
        ukf_var(i,k) = P_UKF(i+8,i+8,k);
    end
end
save('UKF_variance1.mat', 'ukf_var');
%% 
%% === Plot UKF Parameter Estimates (States 9–14) ===
param_labels = {'Tp', 'Kpg', 'Kig', 'Trf', 'Ddn', 'Dup'};
true_params = x_true(9:14, :);
stretch_factor=10;
t_long= t*stretch_factor;
figure('Name','UKF Parameter Estimates','NumberTitle','off');
for i = 1:6
    subplot(3, 2, i);
    plot(t_long(1:Nt-1), xhat_UKF(i+8, 1:Nt-1), 'b-', 'LineWidth', 1.5); hold on;
    plot(t_long(1:Nt-1), true_params(i, 1:Nt-1), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel(param_labels{i});
    title(['Estimate of ', param_labels{i}]);
    legend('Estimated', 'True');
    grid on;
end
sgtitle('UKF Parameter Estimates (Tp, Kpg, Kig, Trf, Ddn, Dup)');
% %% %% === Plot UKF Dynamic State Estimates (States 1–8) ===
% state_labels = {'Vt_filt (x1)', 'PgenFilt (x2)', 'IqPFC (x3)', 'Iq (x4)', ...
%                 'FreqFilt (x5)', 'PI Integral (x6)', 'Pgen (x7)', 'Id (x8)'};
% 
% figure('Name','UKF Dynamic State Estimates','NumberTitle','off');
% for i = 1:8
%     subplot(4, 2, i);
%     plot(t_long(1:Nt-1), xhat_UKF(i, 1:Nt-1), 'b-', 'LineWidth', 1.5); hold on;
%     plot(t_long(1:Nt-1), x_true(i, 1:Nt-1), 'r--', 'LineWidth', 1.5);
%     xlabel('Time (s)'); ylabel(state_labels{i});
%     title(['Estimate vs True: ', state_labels{i}]);
%     legend('Estimated', 'True');
%     grid on;
% end
% sgtitle('UKF Dynamic State Estimates (x1 to x8)');

%% % === SAVE UKF RESULTS ===
UKF_estimates = xhat_UKF;   % UKF output (10 x Nt)

save('UKF_results1.mat', 'UKF_estimates');
disp('✔ UKF results saved to UKF_results1.mat');
