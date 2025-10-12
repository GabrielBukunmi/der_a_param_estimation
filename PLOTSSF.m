
% Parameters
L = -1;
U = 1;
k = 1024;
lambda = (U + L) / 2;
mu = (U - L) / 2;

x = linspace(-2, 2, 1000);
z = (x - lambda) / mu;

z_clip = min(max(z, -1.1), 1.1);

s = z_clip .* (1 + z_clip.^k).^(-1/k);
SSF = lambda + mu .* s;

% Derivative of SSF
dSSF = (1 + z_clip.^k).^(-1 - 1/k);

% Ideal Saturation function
Sat = min(max(x, L), U);

% SDBF and its derivative
SDBF = x - lambda - mu .* s;
dSDBF = dSSF .* ((1 + z_clip.^k).^(1 + 1/k) - 1);

fig1 = figure;

h1 = plot(x, Sat, 'r-', 'LineWidth', 2.0); hold on;

plot(x, SSF, 'b:', 'LineWidth', 1.5);

step = 40;
idx = 2:step:length(x)-1;
idx = [idx length(x)];
plot(x(idx), SSF(idx), 'bo', 'LineWidth', 1.0, 'MarkerSize', 5, ...
     'MarkerFaceColor', 'b', 'HandleVisibility', 'off');

h3 = plot(x, dSSF, '--', 'Color', [0 0.6 0], 'LineWidth', 2.0);

h2 = plot(nan, nan, 'bo:', 'LineWidth', 1.0, 'MarkerSize', 5, ...
          'MarkerFaceColor', 'b');

xlabel('$x$', 'FontSize', 16, 'Interpreter', 'latex');
ylabel('$\mathrm{SSF}\;/\;\mathrm{d}/\mathrm{d}x\,\mathrm{SSF}$', ...
       'FontSize', 16, 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 17);
set(gcf, 'Position', [100 100 850 300]);  % [left bottom width height]

legend([h1, h2, h3], {'Sat', 'SSF $(k=1024)$', '$\mathrm{d}/\mathrm{d}x\,\mathrm{SSF}$'}, ...
       'FontSize', 16, 'Location', 'southeast', 'Interpreter', 'latex');

% ----- PLOT SDBF + DSDBF -----
fig2 = figure;
plot(x, SDBF, 'r-', 'LineWidth', 2); hold on;
plot(x, dSDBF, 'b--', 'LineWidth', 2);
xlabel('$x$', 'FontSize', 17, 'Interpreter', 'latex');
ylabel('$\mathrm{SDBF}\;/\;\mathrm{d}/\mathrm{d}x\,\mathrm{SDBF}$', ...
       'FontSize', 17, 'Interpreter', 'latex');
grid on;
legend('SDBF', '$\mathrm{d}/\mathrm{d}x\,\mathrm{SDBF}$', 'FontSize', 16, ...
       'Location', 'southeast', 'Orientation', 'horizontal', 'Interpreter', 'latex');
set(gca, 'FontSize', 17);
set(gcf, 'Position', [100 100 820 300]);  

% ----- Export SSF Plot -----
exportgraphics(fig1, 'ssf_and_derivative.pdf', 'ContentType', 'vector');
exportgraphics(fig2, 'sdbf_and_derivative.pdf', 'ContentType', 'vector');
