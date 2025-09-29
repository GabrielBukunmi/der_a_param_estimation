% % === Extract DER_A Signals from Simulink Output ===
% Vmag_800   = out.pos_seq_800_Vmag;
% Vmag_802   = out.pos_seq_802_Vmag;
% Imag_800   = out.pos_seq_800_Imag;
% Imag_802   = out.pos_seq_802_Imag;
% Vang_800   = out.pos_seq_800_Vang;
% Vang_802   = out.pos_seq_802_Vang;
% Iang_800   = out.pos_seq_800_Iang;
% Iang_802   = out.pos_seq_802_Iang;
% 
% Vdq0_800   = out.dq0_800_V;
% Vdq0_802   = out.dq0_802_V;
% Idq0_800   = out.dq0_800_I;
% Idq0_802   = out.dq0_802_I;
% 
% f_802      = out.freq_802;
% P_800      = out.pos_seq_800_P;
% Q_800      = out.pos_seq_800_Q;
% Pv_P       = out.pos_seq_PV_P;
% Qv_P       = out.pos_seq_PV_Q;
% IPV_mag =out.pos_seq_I_PV_mag;
% IPV_ang = out.pos_seq_I_PV_ANG;
% IPV_dq0= out.dq0_PV_I_PV;
% 
% % PV2_P      = out.PV850;
% % PV2_Q      = out.Q_PV_850;
% % Id_ref     = out.Id_ref;
% % Iq_ref     =out.Iq_ref;
% % Vd_ref     = out.Vdref;
% % Vq_ref     = out.Vqref;
% 
% 
% % === Save DER_A Signals ===
% save('der_a_DATA2.mat', ...
%     'Vmag_800', 'Vmag_802', ...
%     'Imag_800', 'Imag_802', ...
%     'Vang_800', 'Vang_802', ...
%     'Iang_800', 'Iang_802', ...
%     'Vdq0_800', 'Vdq0_802', ...
%     'Idq0_800', 'Idq0_802', ...
%     'f_802', 'P_800', 'Q_800','Qv_P','Pv_P','IPV_mag','IPV_ang','IPV_dq0');
% 
% disp('✅ Saved DER_A data to der_a_DATA.mat');


% === Extract DER_A Signals from Simulink Output ===
Vmag_800   = out.pos_seq_800_Vmag;
Vmag_802   = out.pos_seq_802_Vmag;
Imag_800   = out.pos_seq_800_Imag;
Imag_802   = out.pos_seq_802_Imag;
Vang_800   = out.pos_seq_800_Vang;
Vang_802   = out.pos_seq_802_Vang;
Iang_800   = out.pos_seq_800_Iang;
Iang_802   = out.pos_seq_802_Iang;

Vdq0_800   = out.dq0_800_V;
Vdq0_802   = out.dq0_802_V;
Idq0_800   = out.dq0_800_I;
Idq0_802   = out.dq0_802_I;

f_802      = out.freq_802;
P_800      = out.pos_seq_800_P;
Q_800      = out.pos_seq_800_Q;
Pv_P       = out.pos_seq_PV_P;
Qv_P       = out.pos_seq_PV_Q;
% IPV_mag = out.pos_seq_I_PV_mag;
% IPV_ang = out.pos_seq_I_PV_ANG;
% IPV_dq0 = out.dq0_PV_I_PV;

P_PV_842 = out.P_PV_842;
P_PV_850 = out.P_PV_850;
P_PV_832 = out.P_PV_832;
P_PV_840 = out.P_PV_840;
P_PV_824 = out.P_PV_824;

Q_PV_824 = out.Q_PV_824;
Q_PV_840 = out.Q_PV_840;
Q_PV_832 = out.Q_PV_832;
Q_PV_842 = out.Q_PV_842;
Q_PV_850 = out.Q_PV_850;

% === Save DER_A Signals ===
save('der_a_DATAnw.mat', ...
    'Vmag_800', 'Vmag_802', ...
    'Imag_800', 'Imag_802', ...
    'Vang_800', 'Vang_802', ...
    'Iang_800', 'Iang_802', ...
    'Vdq0_800', 'Vdq0_802', ...
    'Idq0_800', 'Idq0_802', ...
    'f_802', 'P_800', 'Q_800', ...
    'Pv_P', 'Qv_P', ...
    'P_PV_842','P_PV_850','P_PV_832','P_PV_840','P_PV_824', ...
    'Q_PV_824','Q_PV_840','Q_PV_832','Q_PV_842','Q_PV_850');

disp('✅ Saved DER_A data to der_a_DATAnw.mat');


% % === Extract DER_A Signals from Simulink Output ===
% Vmag_800   = out.pos_seq_800_Vmag;
% Vmag_802   = out.pos_seq_802_Vmag;
% Imag_800   = out.pos_seq_800_Imag;
% Imag_802   = out.pos_seq_802_Imag;
% Vang_800   = out.pos_seq_800_Vang;
% Vang_802   = out.pos_seq_802_Vang;
% Iang_800   = out.pos_seq_800_Iang;
% Iang_802   = out.pos_seq_802_Iang;
% 
% Vdq0_800   = out.dq0_800_V;
% Vdq0_802   = out.dq0_802_V;
% Idq0_800   = out.dq0_800_I;
% Idq0_802   = out.dq0_802_I;
% 
% f_802      = out.freq_802;
% P_800      = out.pos_seq_800_P;
% Q_800      = out.pos_seq_800_Q;
% Pv_P       = out.pos_seq_PV_P;
% Qv_P       = out.pos_seq_PV_Q;
% IPV_mag = out.pos_seq_I_PV_mag;
% IPV_ang = out.pos_seq_I_PV_ANG;
% IPV_dq0 = out.dq0_PV_I_PV;
% 
% Vpv = out.pos_seq_V_PV_mag;  % Added Vpv signal extraction
% 
% % === Save DER_A Signals ===
% save('der_a_DATA_onePV.mat', ...
%     'Vmag_800', 'Vmag_802', ...
%     'Imag_800', 'Imag_802', ...
%     'Vang_800', 'Vang_802', ...
%     'Iang_800', 'Iang_802', ...
%     'Vdq0_800', 'Vdq0_802', ...
%     'Idq0_800', 'Idq0_802', ...
%     'f_802', 'P_800', 'Q_800', ...
%     'Qv_P', 'Pv_P', 'IPV_mag', 'IPV_ang', 'IPV_dq0', ...
%     'Vpv');  % Include Vpv in the saved data
% 
% disp('✅ Saved DER_A data to der_a_DATA2.mat');
