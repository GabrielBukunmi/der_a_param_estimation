clear; clc;

files = {
    struct('path','measurement1_contributions.mat', 'tag','measurement1', 'legendText','$\mathrm{Measurement\ set:}\ \{V, I_{\mathrm{d}}, I_{\mathrm{q}}\}$')
    struct('path','observability_display.mat',       'tag','measurement2', 'legendText','$\mathrm{Measurement\ set:}\ \{V, P, Q\}$')
};

normalize = @(s) regexprep(lower(strtrim(string(s))), '\s+', ' ');
safe = @(s) regexprep(char(s), '[^\w\-]+', '_');
canon_case = @(lbl) local_canon_case(lbl);

for f = 1:numel(files)
    fp  = files{f}.path;
    tag = files{f}.tag;
    legendText = files{f}.legendText;

    S = load(fp,'all_params','all_weights','all_cases');
    params  = S.all_params;
    weights = S.all_weights;
    cases   = S.all_cases;

    cases = string(cases);
    params = string(params);
    norm_cases = normalize(cases);
    labels_norm = unique(norm_cases, 'stable');

    for k = 1:numel(labels_norm)
        sel = (norm_cases == labels_norm(k));
        p = params(sel);
        w = weights(sel);
        if isempty(p), continue; end

        [w, idx] = sort(w, 'descend');
        p = p(idx);
        N = numel(p); y = 0:(N-1);

        ccase = canon_case(labels_norm(k));
        if ccase == "case2"
            lineColor   = [0.2, 0.6, 0.3];
            markerColor = [0, 0, 0];
            lineStyle   = '--';
            markerStyle = 'd';
        else
            lineColor   = [0.2, 0.4, 0.7];
            markerColor = [0.9, 0.2, 0.2];
            lineStyle   = '-';
            markerStyle = 'o';
        end
        lw = 1.4; msz = 52;

        fh = figure('Units','inches','Position',[1 1 7.2 4.2]); hold on;  % fixed size for all
        set(fh,'Renderer','painters');

        for i = 1:N
            line([0, w(i)], [y(i), y(i)], 'Color', lineColor, 'LineWidth', lw, 'LineStyle', lineStyle);
            scatter(w(i), y(i), msz, markerColor, 'filled', 'MarkerEdgeColor', markerColor, 'Marker', markerStyle);
        end

        wmax_all = max(w); if ~isfinite(wmax_all) || wmax_all <= 0, wmax_all = 1; end
        for i = 1:N
            text(w(i) + 0.02*wmax_all, y(i), sprintf('%.2f', w(i)), 'FontSize', 12, ...
                 'HorizontalAlignment','left', 'VerticalAlignment','middle');
        end

        p = cellstr(p);
        p = regexprep(p, '\\text\{([^}]+)\}', '\\mathrm{$1}');
        ax = gca;
        set(ax,'YDir','reverse','TickLabelInterpreter','latex', 'YTick', y, 'YTickLabel', p, 'FontSize', 15);

        if ~contains(tag, 'measurement1')
            xlabel('Weight contributions','FontSize',15);
        end
        ylabel('Parameters','FontSize',15);

        ylim([-0.35, (N-1)+0.35]);
        xmax = max(w); if ~isfinite(xmax) || xmax <= 0, xmax = 1; end
        xlim([0, 1.15*xmax]);

        grid on; box off;

        ax.Units = 'normalized';
        ax.ActivePositionProperty = 'position';
        ax.Position = [0.18 0.14 0.78 0.80];   % same axes box for every figure
        ax.Layer='top'; ax.TickLength=[0 0];
        set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02));

        text(ax, 0.98, 0.04, legendText, ...
            'Units','normalized', 'Interpreter','latex', 'FontSize',13, ...
            'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'Clipping','on');

        outFile = sprintf('%s__%s.pdf', tag, safe(labels_norm(k)));
        set(fh, 'PaperPositionMode', 'auto');
        exportgraphics(fh, outFile, 'ContentType','vector', 'BackgroundColor','none', 'Resolution',600);

        fprintf('Exported: %s\n', outFile);
    end
end

disp('Done.');

function out = local_canon_case(lbl)
    s = lower(strtrim(string(lbl)));
    s = regexprep(s, '\s+', ' ');
    if contains(s, "case 2") || s == "case2"
        out = "case2";
    else
        out = "case1";
    end
end
