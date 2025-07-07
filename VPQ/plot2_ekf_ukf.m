% === LOAD VARIANCE DATA ===
load('EKF_variance.mat', 'ekf_var');
load('UKF_variance.mat', 'ukf_var');
load('EKF_results.mat', 'time_vector');

% === TIME VECTOR ===
Nt = size(ekf_var, 2);
t_plot = time_vector(1:Nt);  % Assumes variance saved for full Nt

% === PARAMETER LABELS ===
param_labels = {'Trv', 'kqv', 'Tg', 'Tiq', 'Tpord'};

% === PLOT ===
figure('Name','EKF vs UKF Parameter Variance','NumberTitle','off');

for i = 1:5
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
