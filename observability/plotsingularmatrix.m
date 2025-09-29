% Plot singular value of observability matrix

data1 = load('obv_data_V_ID_IQ.mat');
data2 = load('obV_data_V_P_Q.mat');

figure;
plot(data1.t, data1.sigma_min, 'b-', 'LineWidth', 2.5); hold on;
plot(data2.t, data2.sigma_min, 'r--', 'LineWidth', 2.5);

xlabel('Time [sec]', 'FontSize', 26);
ylabel('Smallest singular value of $\tilde{O}(t)$', ...
       'Interpreter', 'latex', 'FontSize', 26);

h_legend = legend({'$\{V\ I_d\ I_q\}$', '$\{V\ P\ Q\}$'}, ...
                  'FontSize', 22, 'Interpreter', 'latex', ...
                  'Location', 'best');
title(h_legend, 'Output set', 'FontSize', 22);

grid on;
ax = gca;
ax.FontSize = 20;
ax.Position = [0.15 0.18 0.85 0.73]; 
xlim([0 3.05]);
ylim([0 0.21])
set(gcf, 'Units','inches', 'Position',[1 1 9 5.5]);
   

exportgraphics(gcf, 'observability_comparison.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none', ...
    'Resolution', 600);

%% 

sigma1 = data1.sigma_min(~isnan(data1.sigma_min));
sigma2 = data2.sigma_min(~isnan(data2.sigma_min));

% Compute mean and standard deviation
mean1 = mean(sigma1);
std1  = std(sigma1);

mean2 = mean(sigma2);
std2  = std(sigma2);

% Display results
fprintf('Case: [V Id Iq]  --> Mean = %.4f, STD = %.4f\n', mean1, std1);
fprintf('Case: [V P Q]    --> Mean = %.4f, STD = %.4f\n', mean2, std2);