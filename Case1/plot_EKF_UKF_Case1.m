
addpath(genpath('./Case1_utilities'))
% Load results
load('EKF_results_case1.mat', 'EKF_estimates', 'true_states', 'time_vector');
load('UKF_results_case1.mat', 'UKF_estimates');

t_plot_long = time_vector;

% Match RK4 horizon
Nt = size(true_states, 2);
t_plot_long = t_plot_long(1:Nt-1);

param_labels = {'Trv','kqv','Tg','Tiq','Tpord'};
ylabel_latex = {'$T_{rv}$','$k_{qv}$','$T_{g}$','$T_{iq}$','$T_{pord}$'};

figPos = [100 100 720 405];  
LW      = 3.1;
AXLW    = 0.3;
FSaxis  = 20;
FSlabel = 20;
FSleg   = 16;

for i = 1:5
   
    figPos = [100 100 700 300];  
    
    f = figure('Units','pixels','Position',figPos);
    ax = axes; hold on; grid on; box on
    ax.FontSize  = FSaxis;
    ax.LineWidth = AXLW;

    plot(t_plot_long, true_states(i+5,1:Nt-1), 'k-',  'LineWidth', LW);
    plot(t_plot_long, EKF_estimates(i+5,1:Nt-1),  'b:', 'LineWidth', LW);
    plot(t_plot_long, UKF_estimates(i+5,1:Nt-1),  'r--','LineWidth', LW);

    xlim([t_plot_long(1) t_plot_long(end)]);
    xlabel('Time (s)','Interpreter','latex','FontSize',FSlabel);
    ylabel(ylabel_latex{i},'Interpreter','latex','FontSize',FSlabel);

    legend('Optimal','EKF','UKF','Location','east', ...
           'Interpreter','latex','FontSize',FSleg);

    set(gca,'TickLabelInterpreter','latex');

    fname = sprintf('%s_EKF_UKF.pdf', param_labels{i});
    exportgraphics(gca, fname, 'ContentType','vector', ...
                   'BackgroundColor','white');
end
