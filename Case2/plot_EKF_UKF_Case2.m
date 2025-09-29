
addpath(genpath('./Case2_utilities'))
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

load('EKF_results_case2.mat', 'EKF_estimates', 'true_states', 'time_vector');
load('UKF_results_case2.mat', 'UKF_estimates');

stretch_factor = 1;
t_plot_long = time_vector * stretch_factor;

Nt = size(true_states, 2);
t_plot_long = t_plot_long(1:Nt-1);

param_labels = {'Tp','kpg','Kig','Trf','Ddn','Dup'};
ylabel_disp  = {'$T_p$','$k_{pg}$','$k_{ig}$','$T_{rf}$','Ddn','Dup'};

figPos  = [100 100 720 405];  

LW      = 3.1;
AXLW    = 0.3;
FSaxis  = 20;
FSlabel = 20;
FSleg   = 14;

for i = 1:6
   
    figPos = [100 100 700 300];  
    
    f = figure('Units','pixels','Position',figPos);
    ax = axes; hold on; grid on; box on
    ax.FontSize  = FSaxis;
    ax.LineWidth = AXLW;

    plot(t_plot_long, true_states(i+8,1:Nt-1), 'k-',  'LineWidth', LW);
    plot(t_plot_long, EKF_estimates(i+8,1:Nt-1),  'b:', 'LineWidth', LW);
    plot(t_plot_long, UKF_estimates(i+8,1:Nt-1),  'r--','LineWidth', LW);

    xlim([t_plot_long(1) t_plot_long(end)]);
    xlabel('Time (s)','Interpreter','latex','FontSize',FSlabel);

    if i <= 4
        ylabel(ylabel_disp{i},'Interpreter','latex','FontSize',FSlabel);
    else
        ylabel(ylabel_disp{i},'Interpreter','none','FontSize',FSlabel);
    end

    legend('True','EKF','UKF','Location','east', ...
           'Interpreter','latex','FontSize',FSleg);

    set(gca,'TickLabelInterpreter','latex');

    % export to wide rectangular PDF
    fname = sprintf('%s_EKF_UKF.pdf', param_labels{i});
    exportgraphics(gca, fname, 'ContentType','vector', ...
                   'BackgroundColor','white');
end

