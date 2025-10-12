
font_size = 18;

Y = [0.65, 0.63
     0.61, 0.61
     0.60, 0.21
     0.53, 0.41
     0.48, 0.41
     0.46, 0.38
     0.45, 0.55
     0.45, 0.39
     0.43, 0.42
     0.40, 0.41
     0.19, 0.09
     0.18, 0.08
     0.03, 0.00
     0.01, 0.00
     0.00, 0.00
     0.00, 0.00
     0.00, 0.00];

Ylabel1 = arrayfun(@num2str, Y(:,1), 'UniformOutput', false);
Ylabel1{1} = '';
Ylabel1{2} = '';
Ylabel1{3} = '';
Ylabel1{4} = '';
Ylabel1{5} = '';
Ylabel1{6} = '';
Ylabel1{7} = '';
Ylabel1{8} = '';
Ylabel1{9} = '';
Ylabel1{10} = '';
Ylabel1{15} = '';
Ylabel1{16} = '';
Ylabel1{17} = '';
Ylabel2 = arrayfun(@num2str, Y(:,2), 'UniformOutput', false);
Ylabel2{1} = '';
Ylabel2{2} = '';
Ylabel2{3} = '';
Ylabel2{4} = '';
Ylabel2{5} = '';
Ylabel2{6} = '';
Ylabel2{7} = '';
Ylabel2{8} = '';
Ylabel2{9} = '';
Ylabel2{10} = '';

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
    'YColor',[0 0 0],'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17],'YTickLabel',...
    {'$dP_{\textnormal{max}}$'; '$dP_{\textnormal{min}}$'; '$I_{\textnormal{qmax}}$'; '$I_{\textnormal{qmin}}$'; ...
    '$P_{\textnormal{max}}$'; '$I_{\textnormal{dmax}}$'; '$I_{\textnormal{q$\ell$}}$'; '$P_{\textnormal{min}}$'; ...
    '$I_{\textnormal{q$\mathrm{h}$}}$'; '$I_{\textnormal{dmin}}$'; '$\textnormal{dbd2}$'; '$\textnormal{dbd1}$'; ...
    '$k_{\textnormal{qv}}$'; '$T_{\textnormal{rv}}$'; '$T_{\textnormal{pord}}$'; '$T_{\textnormal{iq}}$'; '$T_{\textnormal{g}}$'});

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'TextColor',[0 0 0],'Interpreter','latex','FontSize',font_size,'EdgeColor',[1 1 1],...
    'Position',[0.55404443764389,0.836247504890724,0.219721226092896,0.07078195861678],...
    'EdgeColor','none','Color','none','Direction','normal');

exportgraphics(figure1,'pl1.pdf','ContentType','vector')