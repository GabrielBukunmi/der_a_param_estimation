clear; clc;

S = load('observability_display.mat','all_params','all_weights','all_cases');
all_params  = S.all_params;
all_weights = S.all_weights;
all_cases   = S.all_cases;

% ---- Sort by descending weight ----
[all_weights, sortIdx] = sort(all_weights,'descend');
all_params  = all_params(sortIdx);
all_cases   = all_cases(sortIdx);

N = numel(all_params);
ypos = 0:(N-1);                 

figure('Units','inches','Position',[1 1 8 max(4.5, 0.28*N)]); hold on;

for i = 1:N
    if strcmp(all_cases{i}, 'case1')
        lineColor   = [0.2, 0.4, 0.7];   markerColor = [0.9, 0.2, 0.2];
        lineStyle   = '-';               markerStyle = 'o';
    else
        lineColor   = [0.2, 0.6, 0.3];   markerColor = [0, 0, 0];
        lineStyle   = '--';              markerStyle = 'd';
    end

    % Plot bar
    line([0, all_weights(i)], [ypos(i), ypos(i)], ...
        'Color', lineColor, 'LineWidth', 1.4, 'LineStyle', lineStyle);

    % Marker
    scatter(all_weights(i), ypos(i), 52, markerColor, 'filled', ...
        'MarkerEdgeColor', markerColor, 'Marker', markerStyle);

    % ---- Add numeric value ----
    text(all_weights(i)+0.02*max(all_weights), ypos(i), ...
        sprintf('%.2f', all_weights(i)), ...
        'FontSize', 12, 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle');
end

set(gca,'YDir','reverse', ...
        'TickLabelInterpreter','latex', ...
        'YTick',ypos, 'YTickLabel',all_params, ...
        'FontSize',15);  
xlabel('Weight contribution','FontSize',15);
ylabel('Parameters','FontSize',15);
ylim([-0.35, (N-1)+0.35]);
xlim([0, max(all_weights)*1.15]);
grid on; box off;
ax = gca; 
ax.Position = [0.15 0.08 0.80 0.88];   
ax.Layer    = 'top';
ax.TickLength = [0 0];

% Legend
case1_dummy = line(NaN,NaN,'LineStyle','-','Color',[0.2,0.4,0.7], ...
                   'Marker','o','MarkerFaceColor',[0.9,0.2,0.2], ...
                   'MarkerEdgeColor',[0.9,0.2,0.2],'LineWidth',1.4,'MarkerSize',6);
case2_dummy = line(NaN,NaN,'LineStyle','--','Color',[0.2,0.6,0.3], ...
                   'Marker','d','MarkerFaceColor',[0,0,0], ...
                   'MarkerEdgeColor',[0,0,0],'LineWidth',1.4,'MarkerSize',6);
legend([case1_dummy, case2_dummy], {'Case 1','Case 2'}, ...
       'Location','southeast','FontSize',13,'Box','off');

set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
exportgraphics(gcf,'weight_contribution.pdf', ...
    'ContentType','vector','BackgroundColor','none','Resolution',600);
