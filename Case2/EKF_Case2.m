%% Pflag =1 Vtrfl=0 Freqflag =1

clear; close all; clc
addpath(genpath('./Case2_utilities'))

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

t = t_unique; 

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
Vbase = 24900/ sqrt(3) * sqrt(2);   
Sbase = 100e3;                      
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


% Power
P800pu_simu = P_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
Q800pu_simu = Q_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
P800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* cos(Vang_800_rad - Iang_800_rad);
Q800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* sin(Vang_800_rad - Iang_800_rad);

% Compute active and reactive power (pu) at node 802
P802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* cos(Vang_802_rad - Iang_802_rad);
Q802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* sin(Vang_802_rad - Iang_802_rad);


%% 

% Sum
P_total = (P_PV_824 + P_PV_842 + P_PV_850 + P_PV_832 + P_PV_840)/100; %Already dividing the data by 1000 in simulink
Q_total = (Q_PV_824 + Q_PV_842 + Q_PV_850 + Q_PV_832 + Q_PV_840)/100;
%% 

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
    'kqv'      ,  6       , ... % proportional voltage control gain (pu/pu)
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
   'Tg'       ,  0.027     , ... % current control time constant
    'Tiq'      ,  0.04     , ... % Q control time constant (s)
    'Tp'       ,  0.018    , ... % transducer time constant (s)
    'Tpord'    ,  0.025     , ... % power order time constant (s)
    'Trf'      ,  0.02     , ... % transducer time constant (seconds) for frequency measurement (> 0)  
    'Trv'      ,  0.019     , ... % transducer time constant (s) for voltage measurement
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
%
n_states = 8;         % Dynamic states: x1 to x8
n_param  = 6;         % Parameters to estimate: x9 to x14
n_total  = n_states + n_param;

dt = 0.02;            
Nt = length(t);       

x_true = zeros(n_total, Nt);

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
Vt_now   = Vmag_802_pu;       
Pref_now = P_total;          
Qref_now = Q_total;          
freq_now = f_802 / 60;       

%% RK4 Propagation of Dynamics
for k = 2:Nt-1
    xk = x_true(:, k-1);

    k1 = compute_k1_P_f_control(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k2 = compute_k2_P_f_control(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k3 = compute_k3_P_f_control(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k4 = compute_k4_P_f_control(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);

    x_true(:, k) = x_true(:, k-1) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end


%% Estimate of state and measurement

n_states = 8;
n_param  = 6;
n_total  = n_states + n_param;
dt  = 0.02;
Nt  = length(t);
nm=6;
xhat_plus = zeros(n_total, Nt);
xhat_minus = zeros(n_total, Nt);
P_plus = zeros(n_total, n_total, Nt);
P_minus = zeros(n_total, n_total, Nt);

% Initial Conditions
x_true(:,1) = [s0_0; s1_0; s2_0; s3_0; s5_0; s6_0; s7_0; s9_0; p.Tp; p.kpg; p.kig; p.Trf; p.Ddn; p.Dup];
xhat_plus(:,1) = [s0_0; s1_0; s2_0; s3_0; s5_0; s6_0; s7_0; s9_0; p.Tp+0.01; p.kpg*0.95; p.kig-0.03; p.Trf-0.0005; p.Ddn-0.5; p.Dup*0.97];

%  Initial Covariance 
P0 = diag([0.01*ones(1,n_states), 0.1*ones(1,n_param)]);
P_plus(:,:,1) = P0;

% Process and Measurement Noise 
Q = (1e-2)^2 * eye(n_total);
R = (1e-4)^2 * eye(nm);
I = eye(n_total);


y_hat = zeros(nm,1);
y_k   = zeros(nm,1);
%% 
for k = 2:Nt-1
    % Previous state estimate
    xk = xhat_plus(:, k-1);

    % RK4 Integration
    k1 = compute_k1_P_f_control(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k2 = compute_k2_P_f_control(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k3 = compute_k3_P_f_control(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k4 = compute_k4_P_f_control(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);

    % Time Update (Prediction Step)
    xhat_minus(:,k) = xk + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
    xhat_minus(1:n_states, k) = min(max(xhat_minus(1:n_states, k), -1e3), 1e3);
    xhat_minus(n_states+1:end, k) = max(xhat_minus(n_states+1:end, k), 1e-6);
    J1 = compute_J1_P_f_control(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    J2 = compute_Jf2_P_f_control(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    J3 = compute_Jf3_P_f_control(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    J4 = compute_Jf4_P_f_control(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    Fk = I + (1/6) * (J1 + 2*J2 + 2*J3 + J4);

    % Covariance Propagation
    P_minus(:,:,k) = Fk * P_plus(:,:,k-1) * Fk' + Q;
    P_minus(:,:,k) = P_minus(:,:,k) + 1e-6 * eye(n_total);  

% === Measurement Prediction ===
x_pred = xhat_minus(:,k);
x1 = x_pred(1);   
x2 = x_pred(2);   
x4 = x_pred(4);   
x5 = x_pred(5);   
x6 = x_pred(6);   
x7 = x_pred(7);   

Tp   = x_pred(9);    
Kpg  = x_pred(10);   
Kig  = x_pred(11);  
Trf  = x_pred(12);   
Ddn  = x_pred(13);   
Dup  = x_pred(14);   

v=1e-6 *randn(6,1);

y_hat = [x1;
         x1 * x4;
         x5;
         x7];

% 
y_k = [x_true(1,k);                      
       x_true(1,k) * x_true(4,k);       
       x_true(5,k);                      
       x_true(7,k)];                    

% 
dy5_hat = (-1/Trf) * x5 + (1/Trf) * freq_now(k);
y_hat = [y_hat; dy5_hat];

dy5_true = (x_true(5,k) - x_true(5,k-1)) / dt;
y_k = [y_k; dy5_true];

% 
Ferr = p.Freq_ref - x5;
uD_dn = Ferr * Ddn;
uD_up = Ferr * Dup;
uD_hat = uD_dn + uD_up;
y_hat = [y_hat; uD_hat];

uD_true = (p.Freq_ref - x_true(5,k)) * x_true(13) + ...
          (p.Freq_ref - x_true(5,k)) * x_true(14);  
y_k = [y_k; uD_true] + v;

% =
Hk = zeros(6, n_total);


Hk(1,1) = 1;


Hk(2,1) = x4;
Hk(2,4) = x1;


Hk(3,5) = 1;


duD_dDdn = Ferr;
duD_dDup = Ferr;
duD_dx5 = -Ddn - Dup;

uD = uD_dn + uD_up;

A = aux_fcn_SSF(Pref_now(k) + uD - x2, p.femax, p.femin, p.k);
dA_dx2 = -aux_fcn_dSSF(Pref_now(k) + uD - x2, p.femax, p.femin, p.k);
dA_dDdn = aux_fcn_dSSF(Pref_now(k) + uD - x2, p.femax, p.femin, p.k) * duD_dDdn;
dA_dDup = aux_fcn_dSSF(Pref_now(k) + uD - x2, p.femax, p.femin, p.k) * duD_dDup;
dA_dx5  = aux_fcn_dSSF(Pref_now(k) + uD - x2, p.femax, p.femin, p.k) * duD_dx5;

B = aux_fcn_SSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k);
dB_dx6 = aux_fcn_dSSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k);
dB_dkpg = aux_fcn_dSSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k) * A;
dB_dx2 = aux_fcn_dSSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k) * Kpg * dA_dx2;
dB_dDdn = aux_fcn_dSSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k) * Kpg * dA_dDdn;
dB_dDup = aux_fcn_dSSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k) * Kpg * dA_dDup;
dB_dx5 = aux_fcn_dSSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k) * Kpg * dA_dx5;
dx6_dKig = dt * A;
dB_dKig = aux_fcn_dSSF(x6 + Kpg * A, p.Pmax, p.Pmin, p.k) * dx6_dKig;

SSF_x7 = aux_fcn_SSF(x7, p.Pmax, p.Pmin, p.k);
dSSF_x7 = aux_fcn_dSSF(x7, p.Pmax, p.Pmin, p.k);
df_dx7 = (1/Tp + p.kw) * dSSF_x7 - p.kw;

Hk(4,2)  = (1/Tp) * dB_dx2 / df_dx7;      % x2
Hk(4,5)  = (1/Tp) * dB_dx5 / df_dx7;      % x5
Hk(4,6)  = (1/Tp) * dB_dx6 / df_dx7;      % x6
Hk(4,7)  = 1;                              % x7
Hk(4,9)  = ((SSF_x7 - B) / Tp^2) / df_dx7; % Tp
Hk(4,10) = (1/Tp) * dB_dkpg / df_dx7;     % Kpg
Hk(4,13) = (1/Tp) * dB_dDdn / df_dx7;     % Ddn
Hk(4,14) = (1/Tp) * dB_dDup / df_dx7;     % Dup
Hk(4,11) = (1/Tp) * dB_dKig / df_dx7;     % Kig
Hk(5,5)  = -1/Trf;
Hk(5,12) = (x5 - freq_now(k)) / Trf^2;
Hk(6,5)  = -Ddn - Dup;  
Hk(6,13) = Ferr;        
Hk(6,14) = Ferr;        

%  EKF Correction Step
innovation = y_k - y_hat;
S = Hk * P_minus(:,:,k) * Hk' + R;
Kk = P_minus(:,:,k) * Hk' / S;

xhat_plus(:,k) = xhat_minus(:,k) + Kk * innovation;
xhat_plus(1:n_states,k) = min(max(xhat_plus(1:n_states,k), -1e6), 1e6);
P_plus(:,:,k) = (I - Kk * Hk) * P_minus(:,:,k);


end


ekf_var = zeros(6, Nt);
for i = 1:6
    for k = 1:Nt
        ekf_var(i,k) = P_plus(i+8, i+8, k);
    end
end

% Optional: save for reuse
save('EKF_variance1.mat', 'ekf_var');


%% EKF Parameter Estimate Plot
true_params = x_true(9:14, :); 
param_names = {'Tp', 'Kpg', 'Kig', 'Trf', 'Ddn', 'Dup'};
stretch_factor=10;
t_long= t*stretch_factor;
figure;
for i = 1:6
    subplot(3, 2, i);
    plot(t_long(1:Nt-1), xhat_plus(i+8, 1:Nt-1), 'b-', 'LineWidth', 1.5); hold on;
    plot(t_long(1:Nt-1), true_params(i, 1:Nt-1), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel(param_names{i});
    title(['Estimate of ', param_names{i}]);
    legend('Estimated', 'True');
    grid on;
end
sgtitle('EKF Parameter Estimates (Tp, Kpg, Kig, Trf, Ddn, Dup)');

%% %% % === SAVE EKF RESULTS AND TRUE STATES ===
EKF_estimates = xhat_plus;   
true_states   = x_true;     
time_vector   = t;           

save('EKF_results_case2.mat', 'EKF_estimates', 'true_states', 'time_vector');
disp('✔ EKF results saved to EKF_results1.mat');

