%% Обработка численного эксперимента

% Автор: Баталов С. А.
% магистратура, курс 1, 2023

% Схема данных численного эксперимента
% x(1), x/L(2), dx(3), x_CO2(4), n_CO2(5), n_AR(6), rho_CO2(7), rho_AR(8),
% p(9), v(10), T(11), T12(12), T3(13), k(14), a(15), M(16), dV(17), trQ(18), 
% vQ12(19), vQ3(20), dQ(21), tdQ(22), sVisc(23), bVisc(24), xxP(25), L(26),
% dx/L(27)

%% Подготовка данных

% Размер матрицы значений
range = [1, 1, 100, 27];

% Путь расположения необходимых файлов
path = "../../QT/BuildMixtureSolverCO2Ar/results";

% Чтение данных теоретического расчета
mp_1 = readmatrix(path + "/" + "mp_dx100_CO2AR_m5_x50_T3_bV0_dV0_p100.txt", ...
    'Range', range);
mp_2 = readmatrix(path + "/" + "mp_dx200_CO2AR_m5_x50_T3_bV0_dV0_p100.txt", ...
    'Range', range);
mp_3 = readmatrix(path + "/" + "mp_dx400_CO2AR_m5_x50_T3_bV0_dV0_p100.txt", ...
    'Range', range);
mp_4 = readmatrix(path + "/" + "mp_dx100_CO2AR_m5_x50_T3_bV0_dV0_p1000.txt", ...
    'Range', range);
mp_5 = readmatrix(path + "/" + "mp_dx100_CO2AR_m3_x50_T3_bV0_dV0_p1000.txt", ...
    'Range', range);
mp_6 = readmatrix(path + "/" + "mp_dx50_CO2AR_m3_x50_T3_bV0_dV0_p1000.txt", ...
    'Range', range);
mp_7 = readmatrix(path + "/" + "mp_dx25_CO2AR_m3_x50_T3_bV0_dV0_p1000.txt", ...
    'Range', range);
mp_8 = readmatrix(path + "/" + "mp_dx50_CO2AR_m2_x50_T3_bV0_dV0_p1000.txt", ...
    'Range', range);
mp_9 = readmatrix(path + "/" + "mp_dx25_CO2AR_m2_x50_T3_bV0_dV0_p2000.txt", ...
    'Range', range);
mp_10 = readmatrix(path + "/" + "mp_dx25_CO2AR_m6_x50_T3_bV0_dV0_p2000.txt", ...
    'Range', range);
mp_11 = readmatrix(path + "/" + "mp_dx25_CO2AR_m2_x50_T3_bV0_dV1_p2000.txt", ...
    'Range', range);
mp_12 = readmatrix(path + "/" + "mp_dx25_CO2AR_m2_x50_T3_bV1_dV0_p2000.txt", ...
    'Range', range);
mp_13 = readmatrix(path + "/" + "mp_dx25_CO2AR_m2_x50_T3_bV1_dV1_p2000.txt", ...
    'Range', range);

% Цветной и Ч-Б варианты
ls = {'-', '--', '-.', ':'};
col = {[1.0000, 0.0000, 0.0000], ...
       [0.4660, 0.6740, 0.1880], ...
       [0.0000, 0.0000, 1.0000], ...
       [0.3010, 0.7450, 0.9330], ...
       [0.0000, 0.0000, 0.0000]};

% Настроечные параметры
l_width = 1.3;
f_name = 'Times New Roman';
f_size = 20;
pos = [0, 0, 550, 500];
pos_ = [0, 0, 1100, 600];

% Масштаб отображения
x1 = 1:60;
x2 = 1:35;
x3 = 1:25;

% Установка стиля и размера шрифта
set(0, 'DefaultAxesFontSize', f_size, 'DefaultAxesFontName', f_name);
set(0, 'DefaultTextFontSize', f_size, 'DefaultTextFontName', f_name);

%% Группа графиков 1

% График распределения температуры
f1 = figure('OuterPosition', pos_);
p = plot(mp_1(x1, 2), mp_1(x1, 11), mp_2(x2, 2), mp_2(x2, 11), mp_3(x2, 2), mp_3(x2, 11), ...
    mp_1(x1, 2), mp_1(x1, 12), mp_2(x2, 2), mp_2(x2, 12), mp_3(x2, 2), mp_3(x2, 12), ...
    mp_1(x1, 2), mp_1(x1, 13), mp_2(x2, 2), mp_2(x2, 13), mp_3(x2, 2), mp_3(x2, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{1}; p(2).Color = col{2};
p(3).LineStyle = ls{1}; p(3).Color = col{3};
p(4).LineStyle = ls{2}; p(4).Color = col{1};
p(5).LineStyle = ls{2}; p(5).Color = col{2};
p(6).LineStyle = ls{2}; p(6).Color = col{3};
p(7).LineStyle = ls{3}; p(7).Color = col{1};
p(8).LineStyle = ls{3}; p(8).Color = col{2};
p(9).LineStyle = ls{3}; p(9).Color = col{3};
title('100 Па, 300 К, M = 5, x = 0.5')
l = legend('$T$', '$T$', '$T$', '$T_{12}$', '$T_{12}$', '$T_{12}$', ...
    '$T_{3},\,dx_0 = 100\,\mu m$', '$T_{3},\,dx_0 = 200\,\mu m$', ...
    '$T_{3},\,dx_0 = 400\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f1, 'figures/colored/rus/fig_T_p100_M5.jpg');

% График распределения длин ячеек в метрах
f2 = figure('OuterPosition', pos);
p = plot(mp_1(x2, 2), mp_1(x2, 3) * 1e6, ...
    mp_2(x2, 2), mp_2(x2, 3) * 1e6, mp_3(x2, 2), mp_3(x2, 3) * 1e6, ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{1}; p(2).Color = col{2};
p(3).LineStyle = ls{1}; p(3).Color = col{3};
title('100 Па, 300 К, M = 5, x = 0.5')
l = legend('$dx_0 = 100\,\mu m$', '$dx_0 = 200\,\mu m$', '$dx_0 = 400\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f2, 'figures/colored/rus/fig_dx_p100_M5.jpg');

% График распределения длин ячеек в длинах свободного пробега
f3 = figure('OuterPosition', pos);
p = plot(mp_1(x2, 2), mp_1(x2, 27), mp_2(x2, 2), mp_2(x2, 27), ...
    mp_3(x2, 2), mp_3(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{1}; p(2).Color = col{2};
p(3).LineStyle = ls{1}; p(3).Color = col{3};
title('100 Па, 300 К, M = 5, x = 0.5')
l = legend('$dx_0 = 100\,\mu m$', '$dx_0 = 200\,\mu m$', '$dx_0 = 400\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f3, 'figures/colored/rus/fig_dx_L_p100_M5.jpg');

%% Группа графиков 2

% График распределения температуры
f4 = figure('OuterPosition', pos_);
p = plot(mp_4(x2, 2), mp_4(x2, 11), mp_4(x2, 2), mp_4(x2, 12), mp_4(x2, 2), mp_4(x2, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{2};
p(3).LineStyle = ls{3}; p(3).Color = col{3};
title('1000 Па, 300 К, M = 5, x = 0.5')
l = legend('$T$', '$T_{12}$', '$T_{3},\,dx_0 = 100\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f4, 'figures/colored/rus/fig_T_p1000_M5_x100.jpg');

% График распределения длин ячеек в метрах
f5 = figure('OuterPosition', pos);
p = plot(mp_4(x2, 2), mp_4(x2, 3) * 1e6, 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 5, x = 0.5')
l = legend('$dx_0 = 100\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f5, 'figures/colored/rus/fig_dx_p1000_M5_x100.jpg');

% График распределения длин ячеек в длинах свободного пробега
f6 = figure('OuterPosition', pos);
p = plot(mp_4(x2, 2), mp_4(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 5, x = 0.5')
l = legend('$dx_0 = 100\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f6, 'figures/colored/rus/fig_dx_L_p1000_M5_x100.jpg');

%% Группа графиков 3

% График распределения температуры
f7 = figure('OuterPosition', pos_);
p = plot(mp_5(x2, 2), mp_5(x2, 11), mp_5(x2, 2), mp_5(x2, 12), mp_5(x2, 2), mp_5(x2, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{2};
p(3).LineStyle = ls{3}; p(3).Color = col{3};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$T$', '$T_{12}$', '$T_{3},\,dx_0 = 100\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f7, 'figures/colored/rus/fig_T_p1000_M3_x100.jpg');

% График распределения длин ячеек в метрах
f8 = figure('OuterPosition', pos);
p = plot(mp_5(x2, 2), mp_5(x2, 3) * 1e6, 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 100\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f8, 'figures/colored/rus/fig_dx_p1000_M3_x100.jpg');

% График распределения длин ячеек в длинах свободного пробега
f9 = figure('OuterPosition', pos);
p = plot(mp_5(x2, 2), mp_5(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 100\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f9, 'figures/colored/rus/fig_dx_L_p1000_M3_x100.jpg');

%% Группа графиков 4

% График распределения температуры
f10 = figure('OuterPosition', pos_);
p = plot(mp_6(x2, 2), mp_6(x2, 11), mp_6(x2, 2), mp_6(x2, 12), mp_6(x2, 2), mp_6(x2, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{2};
p(3).LineStyle = ls{3}; p(3).Color = col{3};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$T$', '$T_{12}$', '$T_{3},\,dx_0 = 50\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f10, 'figures/colored/rus/fig_T_p1000_M3_x50.jpg');

% График распределения длин ячеек в метрах
f11 = figure('OuterPosition', pos);
p = plot(mp_6(x2, 2), mp_6(x2, 3) * 1e6, 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 50\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f11, 'figures/colored/rus/fig_dx_p1000_M3_x50.jpg');

% График распределения длин ячеек в длинах свободного пробега
f12 = figure('OuterPosition', pos);
p = plot(mp_6(x2, 2), mp_6(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 50\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f12, 'figures/colored/rus/fig_dx_L_p1000_M3_x50.jpg');

%% Группа графиков 5

% График распределения температуры
f13 = figure('OuterPosition', pos_);
p = plot(mp_7(x2, 2), mp_7(x2, 11), mp_7(x2, 2), mp_7(x2, 12), mp_7(x2, 2), mp_7(x2, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{2};
p(3).LineStyle = ls{3}; p(3).Color = col{3};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$T$', '$T_{12}$', '$T_{3},\,dx_0 = 25\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f13, 'figures/colored/rus/fig_T_p1000_M3_x25.jpg');

% График распределения длин ячеек в метрах
f14 = figure('OuterPosition', pos);
p = plot(mp_7(x2, 2), mp_7(x2, 3) * 1e6, 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f14, 'figures/colored/rus/fig_dx_p1000_M3_x25.jpg');

% График распределения длин ячеек в длинах свободного пробега
f15 = figure('OuterPosition', pos);
p = plot(mp_7(x2, 2), mp_7(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f15, 'figures/colored/rus/fig_dx_L_p1000_M3_x25.jpg');

%% Группа графиков 6

% График распределения температуры
f16 = figure('OuterPosition', pos_);
p = plot(mp_7(x2, 2), mp_7(x2, 11), mp_6(x2, 2), mp_6(x2, 11), mp_5(x3, 2), mp_5(x3, 11), ...
    mp_7(x2, 2), mp_7(x2, 12), mp_6(x2, 2), mp_6(x2, 12), mp_5(x3, 2), mp_5(x3, 12), ...
    mp_7(x2, 2), mp_7(x2, 13), mp_6(x2, 2), mp_6(x2, 13), mp_5(x3, 2), mp_5(x3, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{1}; p(2).Color = col{2};
p(3).LineStyle = ls{1}; p(3).Color = col{3};
p(4).LineStyle = ls{2}; p(4).Color = col{1};
p(5).LineStyle = ls{2}; p(5).Color = col{2};
p(6).LineStyle = ls{2}; p(6).Color = col{3};
p(7).LineStyle = ls{3}; p(7).Color = col{1};
p(8).LineStyle = ls{3}; p(8).Color = col{2};
p(9).LineStyle = ls{3}; p(9).Color = col{3};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$T$', '$T$', '$T$', '$T_{12}$', '$T_{12}$', '$T_{12}$', ...
    '$T_{3},\,dx_0 = 25\,\mu m$', '$T_{3},\,dx_0 = 50\,\mu m$', ...
    '$T_{3},\,dx_0 = 100\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
%saveas(f16, 'figures/colored/rus/fig_T_p1000_M3.jpg');

% График распределения длин ячеек в метрах
f17 = figure('OuterPosition', pos);
p = plot(mp_7(x2, 2), mp_7(x2, 3) * 1e6, ...
    mp_6(x2, 2), mp_6(x2, 3) * 1e6, mp_5(x2, 2), mp_5(x2, 3) * 1e6, ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{1}; p(2).Color = col{2};
p(3).LineStyle = ls{1}; p(3).Color = col{3};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', '$dx_0 = 50\,\mu m$', '$dx_0 = 100\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f17, 'figures/colored/rus/fig_dx_p1000_M3.jpg');

% График распределения длин ячеек в длинах свободного пробега
f18 = figure('OuterPosition', pos);
p = plot(mp_7(x2, 2), mp_7(x2, 27), mp_6(x2, 2), mp_6(x2, 27), ...
    mp_5(x2, 2), mp_5(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{1}; p(2).Color = col{2};
p(3).LineStyle = ls{1}; p(3).Color = col{3};
title('1000 Па, 300 К, M = 3, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', '$dx_0 = 50\,\mu m$', '$dx_0 = 100\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f18, 'figures/colored/rus/fig_dx_L_p1000_M3.jpg');

%% Группа графиков 7

% График распределения температуры
f19 = figure('OuterPosition', pos_);
p = plot(mp_8(x2, 2), mp_8(x2, 11), mp_8(x2, 2), mp_8(x2, 12), mp_8(x2, 2), mp_8(x2, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{2};
p(3).LineStyle = ls{3}; p(3).Color = col{3};
title('1000 Па, 300 К, M = 2, x = 0.5')
l = legend('$T$', '$T_{12}$', '$T_{3},\,dx_0 = 50\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f19, 'figures/colored/rus/fig_T_p1000_M2_x50.jpg');

% График распределения длин ячеек в метрах
f20 = figure('OuterPosition', pos);
p = plot(mp_8(x2, 2), mp_8(x2, 3) * 1e6, 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 2, x = 0.5')
l = legend('$dx_0 = 50\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f20, 'figures/colored/rus/fig_dx_p1000_M2_x50.jpg');

% График распределения длин ячеек в длинах свободного пробега
f21 = figure('OuterPosition', pos);
p = plot(mp_8(x2, 2), mp_8(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('1000 Па, 300 К, M = 2, x = 0.5')
l = legend('$dx_0 = 50\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f21, 'figures/colored/rus/fig_dx_L_p1000_M2_x50.jpg');

%% Группа графиков 8

% График распределения температуры
f22 = figure('OuterPosition', pos_);
p = plot(mp_9(x2, 2), mp_9(x2, 11), mp_9(x2, 2), mp_9(x2, 12), mp_9(x2, 2), mp_9(x2, 13), ...
    'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{2};
p(3).LineStyle = ls{3}; p(3).Color = col{3};
title('2000 Па, 300 К, M = 2, x = 0.5')
l = legend('$T$', '$T_{12}$', '$T_{3},\,dx_0 = 25\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f22, 'figures/colored/rus/fig_T_p2000_M2_x25.jpg');

% График распределения длин ячеек в метрах
f23 = figure('OuterPosition', pos);
p = plot(mp_9(x2, 2), mp_9(x2, 3) * 1e6, 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('2000 Па, 300 К, M = 2, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f23, 'figures/colored/rus/fig_dx_p2000_M2_x25.jpg');

% График распределения длин ячеек в длинах свободного пробега
f24 = figure('OuterPosition', pos);
p = plot(mp_9(x2, 2), mp_9(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('2000 Па, 300 К, M = 2, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f24, 'figures/colored/rus/fig_dx_L_p2000_M2_x25.jpg');

%% Группа графиков 9

% График распределения температуры
f25 = figure('OuterPosition', pos_);
p = plot(mp_10(x3, 2), mp_10(x3, 11), mp_10(x3, 2), mp_10(x3, 12), mp_10(x3, 2), ...
    mp_10(x3, 13), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{2};
p(3).LineStyle = ls{3}; p(3).Color = col{3};
title('2000 Па, 300 К, M = 6, x = 0.5')
l = legend('$T$', '$T_{12}$', '$T_{3},\,dx_0 = 25\,\mu m$', ...
    'interpreter', 'latex', 'location', 'best');
l.NumColumns = 3;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f25, 'figures/colored/rus/fig_T_p2000_M6_x25.jpg');

% График распределения длин ячеек в метрах
f26 = figure('OuterPosition', pos);
p = plot(mp_10(x2, 2), mp_10(x2, 3) * 1e6, 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('2000 Па, 300 К, M = 6, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx [мкм]', 'FontName', f_name)
axis padded
grid minor
saveas(f26, 'figures/colored/rus/fig_dx_p2000_M6_x25.jpg');

% График распределения длин ячеек в длинах свободного пробега
f27 = figure('OuterPosition', pos);
p = plot(mp_10(x2, 2), mp_10(x2, 27), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
title('2000 Па, 300 К, M = 6, x = 0.5')
l = legend('$dx_0 = 25\,\mu m$', 'interpreter', 'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Длина ячейки dx/L [-]', 'FontName', f_name)
axis padded
grid minor
saveas(f27, 'figures/colored/rus/fig_dx_L_p2000_M6_x25.jpg');

%% Группа графиков 10

% График распределения температуры
f22 = figure('OuterPosition', pos_);
p = plot(mp_9(x3, 2), mp_9(x3, 11), mp_11(x3, 2), mp_11(x3, 11), ...
    mp_12(x3, 2), mp_12(x3, 11), mp_13(x3, 2), mp_13(x3, 11), 'LineWidth', l_width);
p(1).LineStyle = ls{1}; p(1).Color = col{1};
p(2).LineStyle = ls{2}; p(2).Color = col{1};
p(3).LineStyle = ls{1}; p(3).Color = col{2};
p(4).LineStyle = ls{2}; p(4).Color = col{2};
title('2000 Па, 300 К, M = 2, x = 0.5')
l = legend('$\zeta = 0, V_c = 0$', '$\zeta = 0, V_c \neq 0$', ...
    '$\zeta \neq 0, V_c = 0$', '$\zeta \neq 0, V_c \neq 0$', 'interpreter', ...
    'latex', 'location', 'best');
l.NumColumns = 1;
xlabel('$\overline{x}$', 'Interpreter', 'latex')
ylabel('Температура [К]', 'FontName', f_name)
axis padded
grid minor
saveas(f22, 'figures/colored/rus/fig_T_p2000_M2.jpg');