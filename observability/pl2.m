
font_size = 18;

Y = [0.63, 0.61
     0.49, 0.50
     0.39, 0.35
     0.38, 0.40
     0.36, 0.30
     0.19, 0.10
     0.16, 0.12
     0.01, 0.00
     0.01, 0.00
     0.00, 0.00
     0.00, 0.00
     0.00, 0.00
     0.00, 0.00];

Ylabel1 = arrayfun(@num2str, Y(:,1), 'UniformOutput', false);
Ylabel1{1} = '';
Ylabel1{2} = '';
Ylabel1{3} = '';
Ylabel1{4} = '';
Ylabel1{10} = '';
Ylabel1{11} = '';
Ylabel1{12} = '';
Ylabel1{13} = '';
Ylabel2 = arrayfun(@num2str, Y(:,2), 'UniformOutput', false);
Ylabel2{1} = '';
Ylabel2{2} = '';
Ylabel2{3} = '';
Ylabel2{4} = '';

% Create figure
figure1 = figure('Color',[1 1 1],'Theme','light','Position',[100, 100, 400, 700]);

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to bar
bar1 = bar(Y,'Horizontal','on','FontSize',16,'GroupWidth',0.4,'BarWidth',0.6);
set(bar1(1),'DisplayName','$\{V,\,I_{\textnormal{d}},\,I_{\textnormal{q}}\}$','Labels',Ylabel1,'FaceColor',[0 0 1],'EdgeColor','none');
set(bar1(2),'DisplayName','$\{V,\,P,\,Q\}$','Labels',Ylabel2,'LineWidth',2,'LineStyle',':','FaceColor','none','EdgeColor',[1 0 0]);

box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',font_size,'TickLabelInterpreter','latex','XColor',[0 0 0],...
    'YColor',[0 0 0],'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13],'YTickLabel',...
    {'$k_{\textnormal{w}}$'; '$P_{\textnormal{min}}$'; '$f_{\textnormal{emin}}$'; '$P_{\textnormal{max}}$'; ...
    '$f_{\textnormal{emax}}$'; '$\textnormal{fdb1}$'; '$\textnormal{fdb2}$'; '$\textnormal{D}_{\textnormal{dn}}$'; ...
    '$T_{\textnormal{p}}$'; '$k_{\textnormal{ig}}$'; '$\textnormal{D}_{\textnormal{up}}$'; '$T_{\textnormal{rf}}$'; '$k_{\textnormal{pg}}$';});

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'TextColor',[0 0 0],'Interpreter','latex','FontSize',font_size,'EdgeColor',[1 1 1],...
    'Position',[0.55404443764389,0.836247504890724,0.219721226092896,0.07078195861678],...
    'EdgeColor','none','Color','none','Direction','normal');

exportgraphics(figure1,'pl2.pdf','ContentType','vector')