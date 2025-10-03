t = linspace(0, 3, 1000);
V = ones(size(t));
V(t >= 1   & t < 1.1) = 0.80;
V(t >= 1.1 & t < 2.0) = 0.90;

fig = figure('Units','inches','Position',[0 0 6 2],'Color','w');
ax  = axes(fig); hold(ax,'on');

h_main = plot(ax, t, V, 'r-', 'LineWidth', 2);

h_vl0 = yline(ax, 0.80, '--', 'Color', [1 0 1],   'LineWidth', 1.5);
h_vl1 = yline(ax, 0.90, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

h_tvl0 = xline(ax, 1.0, ':', 'Color', [1 0 1],   'LineWidth', 1.5);
h_tvl1 = xline(ax, 1.1, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);

xlabel(ax, 'Time [sec]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel(ax, 'Voltage',    'Interpreter', 'latex', 'FontSize', 16);
xlim(ax, [0 3]); ylim(ax, [0.75 1.02]);
grid(ax, 'on');
set(ax, 'FontSize', 12, 'TickDir','out');

set(ax, 'OuterPosition',[0 0 1 1], 'Position',[0.08 0.22 0.90 0.80]);
pbaspect(ax, [3.2 1 1]);


leg = legend(ax, [h_main, h_vl0, h_vl1, h_tvl0, h_tvl1], ...
    {'Voltage', ...
     '$V_{\ell{0}}\!-\!0.80\,\mathrm{pu}$', ...
     '$V_{\ell{1}}\!-\!0.90\,\mathrm{pu}$', ...
     '$t_{\ell{0}}\!-\!0.1\,\mathrm{s}$', ...
     '$t_{\ell{1}}\!-\!0.9\,\mathrm{s}$'}, ...
    'FontSize', 14, 'Box', 'on', 'Interpreter','latex', 'Location','east');
set(leg, 'ItemTokenSize', [9 5]);

exportgraphics(fig, 'voltage_waveform_plot.pdf', 'ContentType','vector', 'BackgroundColor','white');
