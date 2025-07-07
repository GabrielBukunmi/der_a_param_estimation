% === LOAD EKF AND UKF RESULTS ===
load('EKF_results.mat', 'EKF_estimates', 'true_states', 'time_vector');
load('UKF_results.mat', 'UKF_estimates');

% === TIME STRETCHING ===
stretch_factor = 10;
t_plot_long = time_vector * stretch_factor;  % Apply stretching

% === TIME SETUP ===
Nt = size(true_states, 2);
t_plot_long = t_plot_long(1:Nt-1);  % Match RK4 horizon

% === PARAMETER LABELS ===
param_labels = {'Trv', 'kqv', 'Tg', 'Tiq', 'Tpord'};

% === PLOT ALL 5 PARAMETERS ===
figure('Name','True vs EKF vs UKF - Parameter Estimates (Stretched Time)','NumberTitle','off');

for i = 1:5
    subplot(3,2,i);
    
    % Plot True
    plot(t_plot_long, true_states(i+5,1:Nt-1), 'k-', 'LineWidth', 1.9); hold on;
    
    % Plot EKF
    plot(t_plot_long, EKF_estimates(i+5,1:Nt-1), 'b-', 'LineWidth', 1.9);
    
    % Plot UKF
    plot(t_plot_long, UKF_estimates(i+5,1:Nt-1), 'r--', 'LineWidth', 1.9);
    
    title(['Parameter: ', param_labels{i}]);
    xlabel('Time (s)');
    ylabel(param_labels{i});
    legend('True', 'EKF', 'UKF');
    grid on;
end

sgtitle('Comparison of Parameter Estimates: True vs EKF vs UKF (Stretched Time)');
