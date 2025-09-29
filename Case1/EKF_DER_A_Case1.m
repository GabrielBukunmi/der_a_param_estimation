%% Pflag =0 Vtrfl=0 Freqflag =0

clear; close all; clc
addpath(genpath('../util'))

%% 

 load('der_a_DATA.mat')
t = Vmag_800.Time;

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

start_idx = find(t >= 0.001, 1);

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

% Power
P800pu_simu = P_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
Q800pu_simu = Q_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
P800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* cos(Vang_800_rad - Iang_800_rad);
Q800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* sin(Vang_800_rad - Iang_800_rad);

% Compute active and reactive power (pu) at node 802
P802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* cos(Vang_802_rad - Iang_802_rad);
Q802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* sin(Vang_802_rad - Iang_802_rad);


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
    'fdbd1'    , -0.0006 , ... % lower frequency control deadband ≤ 0 (pu)
    'fdbd2'    ,  0.0006 , ... % upper frequency control deadband ≥ 0 (pu)
    'femax'    ,  99       , ... % frequency control maximum error ≥ 0 (pu)
    'femin'    , -99       , ... % frequency control minimum error ≤ 0 (pu)
    'Fh'       ,  60.07    , ... % frequency break-point for high frequency cut-out of inverters
    'Fl'       ,  59.93    , ... % frequency break-point for low frequency cut-out of inverters
    'Freq'     ,  1        , ... % frequency
    'Freq_ref' ,  1        , ... % set-point frequency reference
    'inf'      ,  99    , ... % infinity (algebraic variable)
    'neg_inf'  , -999   ,.... % negative infinity (algebraic variable)
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
    'kw'       ,  10      , ... % time constant for anti-windup (algbraic)
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
    'Ts'       ,  0.04     , ... % Evaluation time of input signal (created by the author of the paper for simplification purposes)
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

s0_0 = Vpu; 
% Compute initial iqpfc value (x3)
s1_0 =0.3598;
s1_0_init = s1_0; 

% Compute initial Iqv and current control state (x4)
Iqv_0 = aux_fcn_Iqv(p, s0_0);                      % From deadband-based function
Iq_ref_0 = aux_fcn_SSF(Iqv_0, p.Iqh1, p.Iql1, p.k); % Saturated reference
s2_0 = aux_fcn_SSF(s1_0 - Iq_ref_0, p.Iqmax, p.Iqmin, p.k);  
s3_0 = P_total(1);  
x_Idpre = [s0_0; s3_0];  
Idpre   = aux_fcn_Idpre_vra(p, x_Idpre);
s4_0     = 0.5;

%% 

n_states = 5;
n_param  = 5;
n_total  = n_states + n_param;

dt  = 0.02;
Nt  = length(t);
x_true = zeros(n_total, Nt);

x_true(:,1) = [ s0_0;       % x1: filtered voltage
                s1_0;       % x2: iqpfc
                s2_0;       % x3: s2
                s3_0;       % x4: Pref controller state (pord)
                s4_0;       % x5: Id control state
                p.Trv;      % x6: Voltage filter time constant
                p.kqv;      % x7: Volt-var gain
                p.Tg;       % x8: Iq control time constant
                p.Tiq;      % x9: Reactive power integrator time constant
                p.Tpord     % x10: Power order time constant
              ];


Vt_now   = Vmag_802_pu;
Pref_now = P_total;
Qref_now = Q_total;
freq_now = f_802 / 60;
%% RK4
for k = 2:Nt-1

    % Extract current augmented state
    xk = x_true(:, k-1);

    k1 = compute_k1_vradynamics(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k2 = compute_k2_vradynamics(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k3 = compute_k3_vradynamics(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k4 = compute_k4_vradynamics(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);

   
    x_true(:, k) = x_true(:, k-1) + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% Time vector
tt = t(1:Nt-1);  

state_labels = {
    'x1: V_filt';
    'x2: iqpfc';
    'x3: Iq';
    'x4: Pgen';
    'x5: Id';
    'x6: Trv';
    'x7: kqv';
    'x8: Tg';
    'x9: Tiq';
    'x10: Tpord'
};


%% Estimate of state and measurement
n_states = 5;
n_param  = 5;
n_total  = n_states + n_param;
dt  = 0.02;
Nt  = length(t);


xhat_plus = zeros(n_total, Nt);
xhat_minus = zeros(n_total, Nt);
P_plus = zeros(n_total, n_total, Nt);
P_minus = zeros(n_total, n_total, Nt);

% --- Initial Conditions ---
x_true(:,1) = [s0_0; s1_0; s2_0; s3_0; s4_0; p.Trv; p.kqv; p.Tg; p.Tiq; p.Tpord];
xhat_plus(:,1) = [s0_0; s1_0; s2_0; s3_0; s4_0; 0.01;p.kqv*0.8; p.Tg*0.6; p.Tiq+0.05; p.Tpord+0.02];

% --- Initial Covariance ---
P0 = diag([0.01*ones(1,5), 0.1*ones(1,5)]);
P_plus(:,:,1) = P0;

% --- Process and Measurement Noise ---
Q = (1e-2)^2 * eye(n_total);           % Process noise
R = (1e-4)^2 * eye(3);              
I = eye(n_total);

y_hat = zeros(3,1);
y_k   = zeros(3,1);

%% EKF Loop
%% === EKF Estimation Loop: Case2
for k = 2:Nt-1
  
    xk = xhat_plus(:, k-1);
   

    % --- RK4 State Prediction ---
    k1 = compute_k1_vradynamics(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k2 = compute_k2_vradynamics(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k3 = compute_k3_vradynamics(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    k4 = compute_k4_vradynamics(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);

    xhat_minus(:,k) = xk + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

    
    xhat_minus(6:10,k)  = max(xhat_minus(6:10,k), 1e-6);  

    % --- Jacobian Fk via RK4 Approximation ---
    J1 = compute_J1_vradynamics(xk, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    J2 = compute_Jf2_vradynamics(xk, k1, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    J3 = compute_Jf3_vradynamics(xk, k2, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    J4 = compute_Jf4_vradynamics(xk, k3, k, dt, p, Vt_now, Qref_now, freq_now, Pref_now);
    Fk = I + (1/6)*(J1 + 2*J2 + 2*J3 + J4);

    % --- Covariance Prediction ---
    P_minus(:,:,k) = Fk * P_plus(:,:,k-1) * Fk' + Q;
    P_minus(:,:,k) = P_minus(:,:,k) + 1e-6 * eye(n_total);  

    % === Measurement update
    x_pred = xhat_minus(:,k);
    x1 = x_pred(1); x3 = x_pred(3); x5 = x_pred(5);
    v=1e-6*randn(3,1);
    y_hat = [x1; x1 * x3; x1 * x5];  % [Vt; Q; P]
    y_k   = [x_true(1,k); x_true(1,k)*x_true(3,k); x_true(1,k)*x_true(5,k)]+v;

    % --- Measurement Jacobian Hk ---
    Hk = zeros(3, n_total);
    Hk(1,1) = 1;             
    Hk(2,1) = x3;  Hk(2,3) = x1;   
    Hk(3,1) = x5;  Hk(3,5) = x1;  

    % === Kalman Update ===
    innovation = y_k - y_hat;
    S = Hk * P_minus(:,:,k) * Hk' + R;
   

    Kk = P_minus(:,:,k) * Hk' / S;
    xhat_plus(:,k) = xhat_minus(:,k) + Kk * innovation;
  
    xhat_plus(6:10,k) = max(xhat_plus(6:10,k), 1e-2);

    % --- Covariance Update ---
    P_plus(:,:,k) = (I - Kk * Hk) * P_minus(:,:,k);
end

% === EKF Estimation-Error Variance Tracking ===
ekf_var = zeros(5, Nt);
for i = 1:5
    for k = 1:Nt
        ekf_var(i,k) = P_plus(i+5, i+5, k);
    end
end


%% Compute F now

%% === EKF Parameter Estimate Plot ===
true_params = x_true(6:10, :); 
param_names = {'Trv', 'kqv', 'Tg', 'Tiq', 'Tpord'};
figure;
for i = 1:5
    subplot(3, 2, i);
    plot(t(1:Nt-1), xhat_plus(i+5, 1:Nt-1), 'b-', 'LineWidth', 1.5); hold on;
    plot(t(1:Nt-1), true_params(i, 1:Nt-1), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel(param_names{i});
    title(['Estimate of ', param_names{i}]);
    legend('Estimated', 'True');
    grid on;
end
sgtitle('EKF Parameter Estimates (Trv, kqv, Tg, Tiq, Tpord)');
disp(xhat_plus(6:10, end-1)')


%% % === SAVE EKF RESULTS 
EKF_estimates = xhat_plus;   
true_states   = x_true;      
time_vector   = t;   

save('EKF_results.mat', 'EKF_estimates', 'true_states', 'time_vector');
disp('✔ EKF results saved to EKF_results.mat');

