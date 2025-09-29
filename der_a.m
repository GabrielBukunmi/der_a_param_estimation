
addpath(genpath('./util'))

%% 

 load('der_a_DATA.mat')
t = Vmag_800.Time;

% Remove duplicate timestamps ===
[t_unique, unique_idx] = unique(t, 'stable');

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

% % Power plots
% figure 
% subplot(1,2,1); plot(t, P800pu_simu,'k-'); hold on; plot(t, P800pu_dq0,'r--'); legend('P from Simulink pos. seq. block (pu)','P calculated in MATLAB pos. seq. (pu)')
% subplot(1,2,2); plot(t, Q800pu_simu,'k-'); hold on; plot(t, Q800pu_dq0,'r--'); legend('Q from Simulink pos. seq. block (pu)','Q calculated in MATLAB pos. seq. (pu)')
% sgtitle('Comparison of Simulink vs. MATLAB Positive Sequence Power (P & Q) at Bus 800')

%% 

% Sum
P_total = (P_PV_824 + P_PV_842 + P_PV_850 + P_PV_832 + P_PV_840)/100; %Already dividing the data by 1000 in simulink
Q_total = (Q_PV_824 + Q_PV_842 + Q_PV_850 + Q_PV_832 + Q_PV_840)/100;
Q_total(1)



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
    'neg_inf'  , -99   ,.... % negative infinity (algebraic variable)
    'Idmax'    ,  9.9      , ... % maximum limit of total active current (or Ipmax)
    'Idmin'    , -9.9      , ... % minimum limit of total active current (or Ipmin)
    'Imax'     ,  1.2      , ... % maximum converter current (pu)
    'Iqh1'     ,  10.2        , ... % maximum limit of reactive current injection, p.u.
    'Iql1'     , -10.2        , ... % minimum limit of reactive current injection, p.u.    
    'Iqmax'    ,  4.0      , ... % maximum limit of total reactive current
    'Iqmin'    , -4.0     , ... % minimum limit of total reactive current
    'k'        ,  1024     , ... % smoothing factor
    'kig'      ,  10       , ... % active power control integral gain
    'kpg'      ,  0.5      , ... % active power control proportional gain
    'kqv'      ,  6       , ... % proportional voltage control gain (pu/pu)
    'kw'       ,  10      , ... % time constant for anti-windup (algbraic)
    'Pmax'     ,   7       , ... % maximum power (pu)
    'Pmin'     , 0        , ... % minimum power (pu)
    'Pref'     ,  0.9      , ... % active power reference at point of fault
    'Qref'     ,  0.1   , ... % reactive power reference
    'rrpwr'    ,  0.05    , ... % Power rise ramp rate following a fault > 0 (pu/s)
    'rrdn'     ,  7        , ... % (created by the paper) - Power rise ramp rate down side, if s9>=0, rrdn = -∞, otherwise, rrdn= -rrpwr
    'rrup'     ,  7       , ... % (created by the paper) - Power rise ramp rate up side, if s9>=0, rrdn = rrpw, otherwise, rrup= ∞
    'Tfh'      ,  4.5      , ... % timer for fh
    'Tfl'      ,  4.6      , ... % timer for fl (Tfl > Trf)
    'Tg'       ,  0.06     , ... % current control time constant
    'Tiq'      ,  0.04    , ... % Q control time constant (s)
    'Tp'       ,  0.02    , ... % transducer time constant (s)
    'Tpord'    ,  0.02     , ... % power order time constant (s)
    'Trf'      ,  0.02     , ... % transducer time constant (seconds) for frequency measurement (> 0)  
    'Trv'      ,  0.02     , ... % transducer time constant (s) for voltage measurement
    'Ts'       ,  0.02     , ... % Evaluation time of input signal (created by the author of the paper for simplification purposes)
    'Tv'       ,  0.02     , ... % time constant on the output of the voltage/frequency cut-out
    'tvh0'     ,  0.16     , ... % timer for vh0 point
    'tvh1'     ,  0.16        , ... % timer for vh1 point
    'tvl0'     ,  0.16     , ... % timer for vl0 point
    'tvl1'     ,  2        , ... % timer for vl1 point
    'vh0'      ,  1.2     , ... % voltage break-point for high voltage cut-out of inverters
    'vh1'      ,  0.98     , ... % voltage break-point for high voltage cut-out of inverters
    'vl0'      ,  0.8      , ... % voltage break-point for low voltage cut-out of inverters
    'vl1'      ,  0.91     , ... % voltage break-point for low voltage cut-out of inverters
    'Vpr'      ,  0.7      , ... % voltage below which frequency tripping is disabled
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
x(1) = s0_0;  
x(9) = s8_0;
s9_0 = solveForS9_0(p,x);


%% initial conditions
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

F = @(t, x, k) dera_indexed(t, x, p, Vmag_802_pu, P_total, Q_total, f_802, k);

% Run fixed-step solver
sol = fixed_step_solver(F, x0, t);
sol.Time = sol.t;
sol.Solution = sol.X;

%
scale = 1/10;
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
% === Save simulated results for this case ===
save('DERA_NERC_results.mat', 't', 'Psim_series', 'Qsim_series');


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
%% 

%
Snerc = load('DERA_NERC_results.mat');
Scal  = load('DERA_CAL_results.mat');


P_meas = P_total*scale;      Q_meas = Q_total*scale;
P_NERC = Snerc.Psim_series*scale;  Q_NERC = Snerc.Qsim_series*scale;
P_CAL  = Scal.Psim_series*scale;   Q_CAL  = Scal.Qsim_series*scale;

nMarkers = 45;
idxP = round(linspace(1, numel(Snerc.t), nMarkers));
idxQ = round(linspace(1, numel(Snerc.t), nMarkers));

fig = figure('Name','Fig8 - Power (pu/10)');

subplot(2,1,1); hold on; grid on
plot(t,       P_meas, 'k-',  'LineWidth', 1.7);                          
plot(Snerc.t, P_NERC, 'b--o','LineWidth', 1.6, ...                       
    'MarkerIndices', idxP, 'MarkerSize', 2, 'MarkerFaceColor','b');
plot(Scal.t,  P_CAL,  'r-.', 'LineWidth', 1.6);                        
ylabel('Active Power (pu)');
legend('Test system', 'NERC values', 'Estimated values', 'Location', 'Best');

subplot(2,1,2); hold on; grid on
plot(t,       Q_meas, 'k-',  'LineWidth', 1.7);
plot(Snerc.t, Q_NERC, 'b--o','LineWidth', 1.6, ...
    'MarkerIndices', idxQ, 'MarkerSize', 2, 'MarkerFaceColor','b');
plot(Scal.t,  Q_CAL,  'r-.', 'LineWidth', 1.6);
xlabel('Time (s)'); ylabel('Reactive Power (pu)');
sgtitle('Comparison of Measured, NERC, and Estimated Active/Reactive Power')


exportgraphics(fig, 'power_comparison_combined.pdf', 'ContentType','vector');
%% 
