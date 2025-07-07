%% Pflag =0 Vtrfl=0 Freqflag =0

clear; close all; clc
addpath(genpath('./util'))

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
start_idx = find(t >= 0.2, 1);

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

% positive sequence plots
figure % plot 1
subplot(2,2,1); plot(t, Vmag_800); hold on; plot(t, Vmag_802); legend('Vmag 800 (kV)','Vmag 802 (kV)')
subplot(2,2,2); plot(t, Vang_800); hold on; plot(t, Vang_802); legend('Vang 800 (deg)','Vang 802 (deg)')
subplot(2,2,3); plot(t, Vmag_800 - Vmag_802); legend('Subtraction: Vmag 800 - Vmag 802 (kV)')
subplot(2,2,4); plot(t, Vang_800 - Vang_802); legend('Subtraction: Vang 800 - Vang 802 (deg)')

figure % plot 2
subplot(2,2,1); plot(t, Vmag_800_pu);  hold on; plot(t, Vmag_802_pu);  legend('Vmag 800 (pu)','Vmag 802 (pu)')
subplot(2,2,2); plot(t, Vang_800_rad); hold on; plot(t, Vang_802_rad); legend('Vang 800 (rad)','Vang 802 (rad)')
subplot(2,2,3); plot(t, Vmag_800_pu  - Vmag_802_pu);  legend('Subtraction: Vmag 800 - Vmag 802 (pu)')
subplot(2,2,4); plot(t, Vang_800_rad - Vang_802_rad); legend('Subtraction: Vang 800 - Vang 802 (rad)')

figure % plot 3
subplot(2,2,1); plot(t, Imag_800); hold on; plot(t, Imag_802); legend('Imag 800 (A)','Imag 802 (A)')
subplot(2,2,2); plot(t, Iang_800); hold on; plot(t, Iang_802); legend('Iang 800 (deg)','Iang 802 (deg)')
subplot(2,2,3); plot(t, Imag_800 - Imag_802); legend('Subtraction: Imag 800 - Imag 802 (A)')
subplot(2,2,4); plot(t, Iang_800 - Iang_802); legend('Subtraction: Iang 800 - Iang 802 (deg)')

figure % plot 4
subplot(2,2,1); plot(t, Imag_800_pu);  hold on; plot(t, Imag_802_pu);  legend('Imag 800 (pu)','Imag 802 (pu)')
subplot(2,2,2); plot(t, Iang_800_rad); hold on; plot(t, Iang_802_rad); legend('Iang 800 (rad)','Iang 802 (rad)')
subplot(2,2,3); plot(t, Imag_800_pu  - Imag_802_pu);  legend('Subtraction: Imag 800 - Imag 802 (pu)')
subplot(2,2,4); plot(t, Iang_800_rad - Iang_802_rad); legend('Subtraction: Iang 800 - Iang 802 (rad)')

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

% dq0 plots
figure % plot 5
subplot(2,3,1); plot(t, Vd_800_pu); legend('Vd 800 (pu)')
subplot(2,3,2); plot(t, Vq_800_pu); legend('Vq 800 (pu)')
subplot(2,3,3); plot(t, V0_800_pu); legend('V0 800 (pu)')
subplot(2,3,4); plot(t, Vd_802_pu); legend('Vd 802 (pu)')
subplot(2,3,5); plot(t, Vq_802_pu); legend('Vq 802 (pu)')
subplot(2,3,6); plot(t, V0_802_pu); legend('V0 802 (pu)')

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
figure % plot 6
subplot(2,3,1); plot(t, V_800_mag); hold on; plot(t, Vmag_800_pu); legend('Vmag800 dq0 (pu)','Vmag800 pos. seq. (pu)')
subplot(2,3,2); plot(t, V_800_mag - Vmag_800_pu); legend('Subtraction: Vmag800 dq0 - Vmag800 pos. seq. (pu)')
subplot(2,3,3); plot(t, V_802_mag - Vmag_802_pu); legend('Subtraction: Vmag802 dq0 - Vmag802 pos. seq. (pu)')
subplot(2,3,4); plot(t, I_800_mag); hold on; plot(t, Imag_800_pu); legend('Imag800 dq0 (pu)','Imag800 pos. seq. (pu)')
subplot(2,3,5); plot(t, I_800_mag - Imag_800_pu); legend('Subtraction: Imag800 dq0 - Imag800 pos. seq. (pu)')
subplot(2,3,6); plot(t, I_802_mag - Imag_802_pu); legend('Subtraction: Imag802 dq0 - Imag802 pos. seq. (pu)')

figure % plot 7
subplot(2,3,1); plot(t, V_800_ang); hold on; plot(t, Vang_800_rad); legend('Vang800 dq0 (rad)','Vang800 pos. seq. (rad)')
subplot(2,3,2); plot(t, V_800_ang - Vang_800_rad); legend('Subtraction: Vang800 dq0 - Vang800 pos. seq. (rad)')
subplot(2,3,3); plot(t, V_802_ang - Vang_802_rad); legend('Subtraction: Vang802 dq0 - Vang802 pos. seq. (rad)')
subplot(2,3,4); plot(t, I_800_ang); hold on; plot(t, Iang_800_rad); legend('Iang800 dq0 (rad)','Iang800 pos. seq. (rad)')
subplot(2,3,5); plot(t, I_800_ang - Iang_800_rad); legend('Subtraction: Iang800 dq0 - Iang800 pos. seq. (rad)')
subplot(2,3,6); plot(t, I_802_ang - Iang_802_rad); legend('Subtraction: Iang802 dq0 - Iang802 pos. seq. (rad)')

% Power
P800pu_simu = P_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
Q800pu_simu = Q_800 ./ Sbase; % from Simulink P/Q measurement block (pos. seq.)
P800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* cos(Vang_800_rad - Iang_800_rad);
Q800pu_dq0  = (3/2) .* Vmag_800_pu .* Imag_800_pu .* sin(Vang_800_rad - Iang_800_rad);

% Compute active and reactive power (pu) at node 802
P802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* cos(Vang_802_rad - Iang_802_rad);
Q802pu_dq0 = (3/2) .* Vmag_802_pu .* Imag_802_pu .* sin(Vang_802_rad - Iang_802_rad);


% Power plots
figure % plot 8
subplot(1,2,1); plot(t, P800pu_simu,'k-'); hold on; plot(t, P800pu_dq0,'r--'); legend('P from Simulink pos. seq. block (pu)','P calculated in MATLAB pos. seq. (pu)')
subplot(1,2,2); plot(t, Q800pu_simu,'k-'); hold on; plot(t, Q800pu_dq0,'r--'); legend('Q from Simulink pos. seq. block (pu)','Q calculated in MATLAB pos. seq. (pu)')

%% 

% Sum
P_total = (P_PV_824 + P_PV_842 + P_PV_850 + P_PV_832 + P_PV_840)/100; %Already dividing the data by 1000 in simulink
Q_total = (Q_PV_824 + Q_PV_842 + Q_PV_850 + Q_PV_832 + Q_PV_840)/100;
Q_total(1)



figure % plot 3 – total PV active and reactive power
subplot(2,1,1); plot(t, P_total, 'LineWidth', 1.2); title('Total PV Active Power'); ylabel('P_{total} (kW)'); grid on
subplot(2,1,2); plot(t, Q_total, 'LineWidth', 1.2); title('Total PV Reactive Power'); ylabel('Q_{total} (kVAR)'); grid on
xlabel('Time (s)');
%% 


% Define DER-A parameters 
p = struct(...                   % Reference: 2023 - NERC - Parameterization of the DER_A Model for Aggregate DER.pdf
    'dbd1'     , -0.5     , ... % lower voltage deadband ≤ 0 (pu)
    'dbd2'     ,  0.5     , ... % upper voltage deadband ≥ 0 (pu)
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
    'neg_inf'  , -99   ,.... % negative infinity (algebraic variable)
    'Idmax'    ,  11.9      , ... % maximum limit of total active current (or Ipmax)
    'Idmin'    , -11.9      , ... % minimum limit of total active current (or Ipmin)
    'Imax'     ,  1.2      , ... % maximum converter current (pu)
    'Iqh1'     ,  10.2        , ... % maximum limit of reactive current injection, p.u.
    'Iql1'     , -10.2        , ... % minimum limit of reactive current injection, p.u.    
    'Iqmax'    ,  6.0      , ... % maximum limit of total reactive current
    'Iqmin'    , -6.0     , ... % minimum limit of total reactive current
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
    'Tg'       ,  0.0079     , ... % current control time constant
    'Tiq'      ,  0.02     , ... % Q control time constant (s)
    'Tp'       ,  0.02    , ... % transducer time constant (s)
    'Tpord'    ,  0.025     , ... % power order time constant (s)
    'Trf'      ,  0.02     , ... % transducer time constant (seconds) for frequency measurement (> 0)  
    'Trv'      ,  0.02     , ... % transducer time constant (s) for voltage measurement
    'Ts'       ,  0.04     , ... % Evaluation time of input signal (created by the author of the paper for simplification purposes)
    'Tv'       ,  0.02     , ... % time constant on the output of the voltage/frequency cut-out
    'tvh0'     ,  0.1     , ... % timer for vh0 point
    'tvh1'     ,  0.2        , ... % timer for vh1 point
    'tvl0'     ,  0.05     , ... % timer for vl0 point
    'tvl1'     ,  2        , ... % timer for vl1 point
    'vh0'      ,  1.1     , ... % voltage break-point for high voltage cut-out of inverters
    'vh1'      ,  1.12     , ... % voltage break-point for high voltage cut-out of inverters
    'vl0'      ,  0.5      , ... % voltage break-point for low voltage cut-out of inverters
    'vl1'      ,  0.88     , ... % voltage break-point for low voltage cut-out of inverters
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



% State initialization
s0_0 = Vpu;
s7_0 = P_total(1);
s8_0 = solveForS8_0(p, s7_0);
s1_0 = aux_fcn_SSF(s8_0, p.Pmax, p.Pmin, p.k);
s2_0 = Q_total(1);
s3_0 = aux_fcn_SSF(s2_0-aux_fcn_SSF(aux_fcn_Iqv(p,s0_0),p.Iqh1,p.Iql1,p.k),p.Iqmax,p.Iqmin,p.k);%iq_est(2)
s4_0 = TripVoltageLogic(t(1), s0_0, p.vl0, p.vl1, p.vh0, p.vh1, p.Vrfrac, p.tvl1, p.tvl0, p.tvh1, p.tvh0);
s5_0 = f_802(1) / Fbase;

% The initialization of s6 requires an iterative numerical method
x(2) = s1_0; 
x(6) = s5_0; 
x(7) = 0.2; % Initial guess for s6_0
x_updated = solveForS6_0(p, x);
s6_0 = x_updated(7);


% Ip estimate from current
% Ip estimate from current
x(1) = s0_0;  
x(9) = s8_0;
s9_0 = solveForS9_0(p,x);


%% Plot the initial conditions
% Gather state variable names and initial values
state_names = {
    's0: Vt (pu)', ...
    's1: Pgen_filt (pu)', ...
    's2: Qref/V (pu)', ...
    's3: Iq (pu)', ...
    's4: Trip Signal', ...
    's5: Frequency (pu)', ...
    's6: PI', ...
    's7: ramp rate (pu)', ...
    's8: Pgen (pu)', ...
    's9: Ip (pu)'
};

initial_values = [
    s0_0;
    s1_0;
    s2_0;
    s3_0;
    s4_0;
    s5_0;
    s6_0;
    s7_0;
    s8_0;
    s9_0;
];


% Initial condition vector (manual)
x0 = [
    s0_0;
    s1_0;
    s2_0;
    s3_0;
    s4_0;
    s5_0;
    s6_0;
    s7_0;
    s8_0;
    s9_0
];

% Define ODE function using index-based data access
F = @(t, x, k) dera_indexed(t, x, p, Vmag_802_pu, P_total, Q_total, f_802, k);

% Run fixed-step solver
sol = fixed_step_solver(F, x0, t);


sol.Time = sol.t;
sol.Solution = sol.X;


pl(sol) % plot solution

% Extract the final values (last row) of the solution matrix
finalValues = sol.Solution(:,end);

% Display final values of each state variable
disp('Final values of the state variables at all flags =0:')
for i = 1:length(finalValues)
    fprintf('state var %2d: %8.4f\n',i,finalValues(i))
end

%

%% GRID INTERFACE SECTION 
% === Extract DER_A outputs ===
Id_series = sol.Solution(10, :);  % Ip = active current (pu)
Iq_series = sol.Solution(4, :);
IqPFC_series = sol.Solution(3, :);% Iq = reactive current (pu)
Vt_series = sol.Solution(1, :); 
Pgen_series = sol.Solution(8, :);% Generated Power(pu)


%Simulation results
Vt_sim = Vt_series;
Psim_series = zeros(size(Id_series));
Qsim_series = zeros(size(Iq_series));



% for k = 1:length(Id_series)
%     Vd = Vt_series(k);  
%     Id = Id_series(k);
%     Iq = Iq_series(k);
%     Pgen =Pgen_series(k);
% 
%     Psim_series(k) = Pgen;
%     Qsim_series(k) =  Vd * Iq;
% end

 %calculate V or alternatively comment 410 -427 and uncomment 400 - 407 
Vt_calculated = zeros(size(Id_series));
for k = 1:length(Id_series)
    Id = Id_series(k);
    Iq = Iq_series(k);  
   Vd = Vt_series(k); 

    % Calculate Ed and Eq
    Eq = 0 + Id*p.Xe;
    Ed = Vt_sim(k)-Iq * p.Xe;

    theta = V_802_ang(k);  

    Vt_calculated(k) = abs((Ed + 1i * Eq) * exp(1i * theta));
    Pgen =Pgen_series(k);
  Psim_series(k) = Ed * Id  + Eq*Iq;
    Qsim_series(k) =  Vt_calculated(k)*Iq;
end

% 


%% Voltage Comparison Plot
figure; hold on;
plot(t, Vt_calculated, 'r-', 'LineWidth', 1.5);
plot(t, Vmag_802_pu, 'k:', 'LineWidth', 2);
legend('Simulated', 'Measured', 'Location', 'Best');
xlabel('Time (s)'); ylabel('Voltage (pu)');
title('Voltage Comparison: Simulated vs Measured'); grid on;


%% Active Power Comparison Plot
figure; hold on;
plot(t, Psim_series, 'r-', 'LineWidth', 1.5);
plot(t, P_total, 'k:', 'LineWidth', 2);
legend('Simulated', 'Measured', 'Location', 'Best');
xlabel('Time (s)'); ylabel('Active Power (pu)');
title('Active Power Comparison: Simulated vs Measured'); grid on;


%% Reactive Power Comparison Plot
figure; hold on;
plot(t, Qsim_series, 'r-', 'LineWidth', 1.5);
plot(t, Q_total, 'k:', 'LineWidth', 2);
legend('Simulated', 'Measured', 'Location', 'Best');
xlabel('Time (s)'); ylabel('Reactive Power (pu)');
title('Reactive Power Comparison: Simulated vs Measured'); grid on;


%% 


function dx = dera_indexed(~, x, p, Vmag, P, Q, f, k)
    dx = zeros(10, 1);

    % Inputs at time step k
    Vt_now    = Vmag(k);
    Pref_now  = P(k);
    Qref_now  = Q(k);
    freq_now  = f(k) / 60;

    % Compute RHS of each state
    dx(1) = -1/p.Trv * x(1) + 1/p.Trv * Vt_now;                                                             % s0
    dx(2) = -1/p.Tp  * x(2) + 1/p.Tp  * aux_fcn_SSF(x(9), p.Pmax, p.Pmin, p.k);                         % s1
    dx(3) = -1/p.Tiq * x(3) + 1/p.Tiq * Qref_now/ aux_fcn_SSF(x(1), p.inf, 0.01, p.k);                 % s2

    dx(4) = -1/p.Tg * x(4) + 1/p.Tg * aux_fcn_SSF(x(3) - aux_fcn_SSF(aux_fcn_Iqv(p, x), p.Iqh1, p.Iql1, p.k), p.Iqmax, p.Iqmin, p.k);       % s3

     dx(5) = -1/p.Tv  * x(5) + 1/p.Tv * TripVoltageLogic(k, x(1), p.vl0, p.vl1, p.vh0, p.vh1, p.Vrfrac, p.tvl1, p.tvl0, p.tvh1, p.tvh0);

    dx(6) = -1/p.Trf * x(6) + freq_now/ p.Trf;                                                           % s5
    dx(7) = p.kig * aux_fcn_A(p, x) + p.kw * (aux_fcn_B(p, x) - p.kpg * aux_fcn_A(p, x) - x(7));        % s6
    dx(8) = (-1/p.Ts * x(8) + 1/p.Ts * Pref_now);                       % s7 modified
    dx(9) = 1/p.Tpord * x(8) - 1/p.Tpord * aux_fcn_SSF(x(9), p.Pmax, p.Pmin, p.k) + p.kw * (aux_fcn_SSF(x(9), p.Pmax, p.Pmin, p.k) - x(9));                                 % s8

    dx(10) = (1/p.Tg * aux_fcn_SSF(aux_fcn_Idpre(p, x), p.Idmax, p.Idmin, p.k) - 1/p.Tg * x(10));                                                                        % s9
end





% function dx_s9 = do_s9_piecewise(p,x)
%     Ip_now = x(10);
%     expr_s9 = (1/p.Tg)*aux_fcn_SSF(aux_fcn_Idpre(p,x), p.Idmax, p.Idmin, p.k) ...
%               - (1/p.Tg)*Ip_now;
% 
%     if Ip_now >= 0
%         upperLimit = p.rrpwr;
%         lowerLimit = p.neg_inf;
%     else
%         upperLimit = p.inf;
%         lowerLimit = -p.rrpwr;
%     end
% 
%     dx_s9 = aux_fcn_SSF(expr_s9, upperLimit, lowerLimit, p.k);
% end

