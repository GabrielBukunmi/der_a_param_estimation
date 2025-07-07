clc;
close all;
% === LOAD VARIANCE DATA ===
load('EKF_variance1.mat', 'ekf_var');
load('UKF_variance1.mat', 'ukf_var');
load('EKF_results1.mat', 'time_vector');

% === TIME VECTOR ===
Nt = size(ekf_var, 2);
t_plot = time_vector(1:Nt);  % Assumes variance saved for full Nt

% === PARAMETER LABELS ===
param_labels = {'Tp', 'kpg', 'Kig', 'Trf', 'Ddn','Dup'};

% === PLOT ===
figure('Name','EKF vs UKF Parameter Variance','NumberTitle','off');

for i = 1:6
    subplot(3,2,i);
    
    semilogy(t_plot, ekf_var(i,:), 'b-', 'LineWidth', 1.5); hold on;
    semilogy(t_plot, ukf_var(i,:), 'r--', 'LineWidth', 1.5);
    
    title(['Variance of Parameter: ', param_labels{i}]);
    xlabel('Time (s)');
    ylabel('Variance (log scale)');
    legend('EKF', 'UKF');
    grid on;
end

sgtitle('Comparison of EKF vs UKF Parameter Estimation Variance');