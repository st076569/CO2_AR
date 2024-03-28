%% Обработка численного эксперимента
% Баталов С. А., курс 4, 2023

%% Подготовка данных

% Чтение данных теоретического расчета
% x(1), x/L(2), dx(3), x_CO2(4), n_CO2(5), n_AR(6), rho_CO2(7), rho_AR(8), 
% p(9), v(10), T(11), T12(12), T3(13), k(14), a(15), M(16), dV(17), trQ(18), 
% vQ12(19), vQ3(20), dQ(21), tdQ(22), sVisc(23), bVisc(24), xxP(25)
mp_T1 = readmatrix('mp_CO2AR_m50_x50_T1_bV1_dV1_p6660.txt', 'Range', [1, 1, 100, 25]);
mp_T2 = readmatrix('mp_CO2AR_m50_x50_T2_bV1_dV1_p6660.txt', 'Range', [1, 1, 100, 25]);
mp_T3 = readmatrix('mp_CO2AR_m50_x50_T3_bV1_dV1_p6660.txt', 'Range', [1, 1, 100, 25]);

% Цветной и Ч-Б варианты
ls = {'-', '--', '-.', ':'};
col = {[1, 0, 0], [0.4660, 0.6740, 0.1880], [0, 0, 1], ...
       [0.3010, 0.7450, 0.9330], [0, 0, 0]};

% Настроечные параметры
l_width = 1.3;
f_name = 'Times New Roman';
f_size = 20;
f_size_ = 20;
pos = [0, 0, 550, 500];
pos_ = [0, 0, 1100, 600];
x1 = 1:60;
x2 = 1:35;

set(0,'DefaultAxesFontSize',f_size,'DefaultAxesFontName',f_name);
%set(0,'DefaultTextFontSize',f_size,'DefaultTextFontName',f_name);

%% Группа графиков 1

% График распределения температуры
f1 = figure('OuterPosition', pos_);
p = plot(mp_T1(x1, 2), mp_T1(x1, 11), ...
         mp_T2(x1, 2), mp_T2(x1, 11), mp_T2(x1, 2), mp_T2(x1, 12), ...
         mp_T3(x1, 2), mp_T3(x1, 11), mp_T3(x1, 2), mp_T3(x1, 12), ...
         mp_T3(x1, 2), mp_T3(x1, 13), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{2};
p(3).Color = col{2};
p(4).LineStyle = ls{1};
p(4).Color = col{3};
p(5).LineStyle = ls{2};
p(5).Color = col{3};
p(6).LineStyle = ls{3};
p(6).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 5, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
l = legend('$T,\,1T$', '$T,\,2T$', '$T_{v},\,2T$', '$T,\,3T$', ...
    '$T_{12},\,3T$', '$T_{3},\,3T$', 'interpreter', 'latex', ...
    'location', 'best', 'FontSize', f_size_);
l.NumColumns = 2;
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f1, 'figures/colored/rus/fig_T_iT.jpg');

% График распределения давления
f2 = figure('OuterPosition', pos);
p = plot(mp_T1(x1, 2), mp_T1(x1, 9), mp_T2(x1, 2), mp_T2(x1, 9), ...
    mp_T3(x1, 2), mp_T3(x1, 9), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{1};
p(3).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 5, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$p,\,1T$', '$p,\,2T$', '$p,\,3T$', 'interpreter', 'latex', ...
    'location', 'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Давление [Па]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f2, 'figures/colored/rus/fig_p_iT.jpg');

% График распределения молярной доли
f3 = figure('OuterPosition', pos);
p = plot(mp_T1(x2, 2), 2 * mp_T1(x2, 4), mp_T2(x2, 2), ...
    2 * mp_T2(x2, 4), mp_T3(x2, 2), 2 * mp_T3(x2, 4), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{1};
p(3).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 5, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$1T$', '$2T$', '$3T$', ...
    'interpreter', 'latex', 'location', 'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('$\widetilde{x}_{_{CO_2}}$', 'FontSize', f_size, 'Interpreter', 'latex')
%ylabel('Молярная доля [%]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f3, 'figures/colored/rus/fig_x_iT.jpg');

% График распределения скорости диффузии
f4 = figure('OuterPosition', pos);
p = plot(mp_T1(x2, 2), mp_T1(x2, 17), mp_T2(x2, 2), mp_T2(x2, 17), ...
    mp_T3(x2, 2), mp_T3(x2, 17), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{1};
p(3).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 5, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$V_{_{CO_2}},\,1T$', '$V_{_{CO_2}},\,2T$', '$V_{_{CO_2}},\,3T$', ...
    'interpreter', 'latex', 'location', 'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Скорость диффузии [м / с]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f4, 'figures/colored/rus/fig_dV_iT.jpg');

% График распределения тепловых потоков
maxQ_T1 = max(abs(mp_T1(x2, 18))) / 100;
maxQ_T2 = max(abs(mp_T2(x2, 18))) / 100;
maxQ_T3 = max(abs(mp_T3(x2, 18))) / 100;
f5 = figure('OuterPosition', pos);
p = plot(mp_T1(1:40, 2), mp_T1(1:40, 21) / maxQ_T1, ...
    mp_T2(1:40, 2), mp_T2(1:40, 21) / maxQ_T2, mp_T3(1:40, 2), ...
    mp_T3(1:40, 21) / maxQ_T3, 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{1};
p(3).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 5, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$\widetilde{q}_{diff},\,1T$', '$\widetilde{q}_{diff},\,2T$', ...
    '$\widetilde{q}_{diff},\,3T$', 'interpreter', 'latex', 'location', ...
    'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('$\widetilde{q}\,\,[\%]$', 'FontSize', f_size, 'Interpreter', 'latex')
%ylabel('Тепловой поток [Дж / (с \cdot м^2)]', 'FontSize', f_size, ...
%    'FontName', f_name)
axis padded
grid minor
saveas(f5, 'figures/colored/rus/fig_Q_iT.jpg');

% График распределения отношения вязкостей
f6 = figure('OuterPosition', pos);
p = plot(mp_T1(x1, 2), mp_T1(x1, 24) ./ mp_T1(x1, 23), mp_T2(x1, 2), ...
    mp_T2(x1, 24) ./ mp_T2(x1, 23), mp_T3(x1, 2), ...
    mp_T3(x1, 24) ./ mp_T3(x1, 23), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{1};
p(3).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 5, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$1T$', '$2T$', '$3T$', 'interpreter', 'latex', 'location', 'best', ...
    'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('\zeta / \eta', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f6, 'figures/colored/rus/fig_bsVisc_iT.jpg');

% График распределения напряжений
f7 = figure('OuterPosition', pos);
p = plot(mp_T1(x2, 2), mp_T1(x2, 25) - mp_T1(x2, 9), ...
    mp_T2(x2, 2), mp_T2(x2, 25) - mp_T2(x2, 9), ...
    mp_T3(x2, 2), mp_T3(x2, 25) - mp_T3(x2, 9), 'LineWidth', l_width);
p(1).LineStyle = ls{1};
p(1).Color = col{1};
p(2).LineStyle = ls{1};
p(2).Color = col{2};
p(3).LineStyle = ls{1};
p(3).Color = col{3};
title('p = 66.6 Па, T = 300 К, M = 5, x_{CO2} = 0.5', 'FontSize', f_size, ...
    'FontWeight', 'normal', 'FontName', f_name)
legend('$\widetilde{P}_{xx},\,1T$', '$\widetilde{P}_{xx},\,2T$', ...
    '$\widetilde{P}_{xx},\,3T$', 'interpreter', 'latex', 'location', ...
    'best', 'FontSize', f_size_)
xlabel('$\overline{x}$', 'FontSize', f_size, 'Interpreter', 'latex')
ylabel('Напряжение [Па]', 'FontSize', f_size, 'FontName', f_name)
axis padded
grid minor
saveas(f7, 'figures/colored/rus/fig_Pxx_iT.jpg')