%% Обработка численного эксперимента
% Баталов С. А., курс 4, 2023

%% Подготовка данных

% Чтение данных теоретического расчета
% x(1), x/L(2), dx(3), x_CO2(4), n_CO2(5), n_AR(6), rho_CO2(7), rho_AR(8), 
% p(9), v(10), T(11), T12(12), T3(13), k(14), a(15), M(16), dV(17), trQ(18), 
% vQ12(19), vQ3(20), dQ(21), tdQ(22), sVisc(23), bVisc(24), Pxx(25)
mp_bV1_dV1 = readmatrix('mp_CO2AR_m70_x50_T3_bV1_dV1_p666.txt', 'Range', [1, 1, 100, 25]);
mp_bV0_dV1 = readmatrix('mp_CO2AR_m70_x50_T3_bV0_dV1_p666.txt', 'Range', [1, 1, 100, 25]);
mp_bV1_dV0 = readmatrix('mp_CO2AR_m70_x50_T3_bV1_dV0_p666.txt', 'Range', [1, 1, 100, 25]);
mp_bV0_dV0 = readmatrix('mp_CO2AR_m70_x50_T3_bV0_dV0_p666.txt', 'Range', [1, 1, 100, 25]);

% Цветной и Ч-Б варианты
ls = {'-', '--', '-.', ':'};
col = {[0.9, 0.2, 0.2], [0.4660, 0.6740, 0.1880], [0, 0, 1], ...
       [0.3010, 0.7450, 0.9330], [0, 0, 0]};

% Настроечные параметры
l_width = 1.3;
f_name = 'Times New Roman';
f_size = 20;
f_size_ = 20;
pos = [0, 0, 550, 500];
pos_ = [0, 0, 1100, 600];
x1 = 1:25;
x2 = 1:35;

set(0,'DefaultAxesFontSize',f_size,'DefaultAxesFontName',f_name);

%% График распределения температуры
f1 = figure('OuterPosition', pos_);
p = plot(mp_bV1_dV1(x1, 2), mp_bV1_dV1(x1, 11), ...
    mp_bV0_dV1(x1, 2), mp_bV0_dV1(x1, 11), ...
    mp_bV1_dV0(x1, 2), mp_bV1_dV0(x1, 11), ...
    mp_bV0_dV0(x1, 2), mp_bV0_dV0(x1, 11), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{3};
p(3).LineStyle = ls{2};
p(3).Color = col{1};
p(4).LineStyle = ls{2};
p(4).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$T,\,\zeta \ne 0,\,V_{_{CO_2}} \ne 0$', ...
    '$T,\,\zeta = 0,\,V_{_{CO_2}} \ne 0$', ...
    '$T,\,\zeta \ne 0,\,V_{_{CO_2}} = 0$', ...
    '$T,\,\zeta = 0,\,V_{_{CO_2}} = 0$', 'interpreter', 'latex', ...
    'location', 'best', 'FontSize', f_size_);
%xlabel('x / L', 'FontSize', f_size, 'FontName', f_name)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f1, 'figures/colored/rus/fig_T_ibVdV.jpg')

%% График распределения давления
f2 = figure('OuterPosition', pos);
p = plot(mp_bV1_dV1(x1, 2), mp_bV1_dV1(x1, 9), ...
    mp_bV0_dV1(x1, 2), mp_bV0_dV1(x1, 9), ...
    mp_bV1_dV0(x1, 2), mp_bV1_dV0(x1, 9), ...
    mp_bV0_dV0(x1, 2), mp_bV0_dV0(x1, 9), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{3};
p(3).LineStyle = ls{2};
p(3).Color = col{1};
p(4).LineStyle = ls{2};
p(4).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$p,\,\zeta \ne 0,\,V_{_{CO_2}} \ne 0$', ...
    '$p,\,\zeta = 0,\,V_{_{CO_2}} \ne 0$', ...
    '$p,\,\zeta \ne 0,\,V_{_{CO_2}} = 0$', ...
    '$p,\,\zeta = 0,\,V_{_{CO_2}} = 0$', 'interpreter', 'latex', ...
    'location', 'best', 'FontSize', f_size_);
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Давление [Па]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f2, 'figures/colored/rus/fig_p_ibVdV.jpg')

%% График распределения молярной доли
f3 = figure('OuterPosition', pos);
p = plot(mp_bV1_dV1(x2, 2), 2 * mp_bV1_dV1(x2, 4), ...
    mp_bV0_dV1(x2, 2), 2 * mp_bV0_dV1(x2, 4), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$\zeta \ne 0$', '$\zeta = 0$', ...
    'interpreter', 'latex', 'location', 'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('$\widetilde{x}_{_{CO_2}}$', 'FontSize', f_size, 'Interpreter', 'latex')
axis padded
grid minor
saveas(f3, 'figures/colored/rus/fig_x_ibVdV.jpg')

%% График распределения скорости диффузии
f4 = figure('OuterPosition', pos);
p = plot(mp_bV1_dV1(x2, 2), mp_bV1_dV1(x2, 17), ...
    mp_bV0_dV1(x2, 2), mp_bV0_dV1(x2, 17), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$V_{_{CO_2}},\,\zeta \ne 0$', '$V_{_{CO_2}},\,\zeta = 0$', ...
    'interpreter', 'latex', 'location', 'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Скорость диффузии [м / с]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f4, 'figures/colored/rus/fig_dV_ibVdV.jpg')

%% График распределения тепловых потоков
maxQ_bV1 = max(abs(mp_bV1_dV1(x2, 18))) / 100;
maxQ_bV0 = max(abs(mp_bV0_dV1(x2, 18))) / 100;
f5 = figure('OuterPosition', pos);
p = plot(mp_bV1_dV1(x2, 2), mp_bV1_dV1(x2, 21) / maxQ_bV1, ...
    mp_bV0_dV1(x2, 2), mp_bV0_dV1(x2, 21) / maxQ_bV0, 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$\zeta \ne 0$', '$\zeta = 0$', 'interpreter', 'latex', ...
    'location', 'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('$\widetilde{q}\,\,[\%]$', 'FontSize', f_size, 'Interpreter', 'latex')
axis padded
grid minor
saveas(f5, 'figures/colored/rus/fig_Q_ibVdV.jpg')

%% График распределения отношения вязкостей
f6 = figure('OuterPosition', pos);
p = plot(mp_bV1_dV1(x1, 2), mp_bV1_dV1(x1, 24) ./ mp_bV1_dV1(x1, 23), ...
    mp_bV1_dV0(x1, 2), mp_bV1_dV0(x1, 24) ./ mp_bV1_dV0(x1, 23), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$V_{_{CO_2}} \ne 0$', '$V_{_{CO_2}} = 0$', ...
    'interpreter', 'latex', 'location', 'best', 'FontSize', f_size_);
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('\zeta / \eta', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f6, 'figures/colored/rus/fig_bsVisc_ibVdV.jpg')

%% График распределения напряжений
f7 = figure('OuterPosition', pos);
p = plot(mp_bV1_dV1(x2, 2), mp_bV1_dV1(x2, 25) - mp_bV1_dV1(x2, 9), ...
    mp_bV0_dV1(x2, 2), mp_bV0_dV1(x2, 25) - mp_bV0_dV1(x2, 9), ...
    mp_bV1_dV0(x2, 2), mp_bV1_dV0(x2, 25) - mp_bV1_dV0(x2, 9), ...
    mp_bV0_dV0(x2, 2), mp_bV0_dV0(x2, 25) - mp_bV0_dV0(x2, 9), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{1};
p(3).Color = col{3};
p(4).LineStyle = ls{1};
p(4).Color = col{4};
title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$\widetilde{P}_{xx},\,\zeta \ne 0,\,V_{_{CO_2}} \ne 0$', ...
    '$\widetilde{P}_{xx},\,\zeta = 0,\,V_{_{CO_2}} \ne 0$', ...
    '$\widetilde{P}_{xx},\,\zeta \ne 0,\,V_{_{CO_2}} = 0$', ...
    '$\widetilde{P}_{xx},\,\zeta = 0,\,V_{_{CO_2}} = 0$', 'interpreter', 'latex', ...
    'location', 'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Напряжение [Па]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f7, 'figures/colored/rus/fig_Pxx_ibVdV.jpg')

% % График распределения показателя адиабаты
% f8 = figure('OuterPosition', pos);
% p = plot(mp_bV1_dV1(x1, 2), mp_bV1_dV1(x1, 14), ...
%     mp_bV0_dV1(x1, 2), mp_bV0_dV1(x1, 14), ...
%     mp_bV1_dV0(x1, 2), mp_bV1_dV0(x1, 14), ...
%     mp_bV0_dV0(x1, 2), mp_bV0_dV0(x1, 14), ...
%     'LineWidth', l_width);
% p(1).LineStyle = ls{1};
% p(1).Color = col{1};
% p(2).LineStyle = ls{1};
% p(2).Color = col{2};
% p(3).LineStyle = ls{1};
% p(3).Color = col{3};
% p(4).LineStyle = ls{1};
% p(4).Color = col{4};
% title('p = 66.6 Па, T = 300 К, M = 7, x_{CO2} = 0.5, 3T', 'FontSize', f_size, ...
%     'FontWeight', 'normal', 'FontName', f_name)
% legend('$k,\,\zeta \ne 0,\,V_{_{CO_2}} \ne 0$', ...
%     '$k,\,\zeta = 0,\,V_{_{CO_2}} \ne 0$', ...
%     '$k,\,\zeta \ne 0,\,V_{_{CO_2}} = 0$', ...
%     '$k,\,\zeta = 0,\,V_{_{CO_2}} = 0$', 'interpreter', 'latex', ...
%     'location', 'best', 'FontSize', f_size_)
% xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
% ylabel('Показатель адиабаты', 'FontSize', f_size, 'FontName', f_name)
% axis padded
% grid minor
% exportgraphics(f8, 'figures/colored/rus/fig_k_ibVdV.eps', 'Resolution', 600)
% exportgraphics(f8, 'figures/colored/rus/fig_k_ibVdV.jpg', 'Resolution', 600)
