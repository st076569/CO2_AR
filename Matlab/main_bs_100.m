%% Обработка численного эксперимента
% Баталов С. А., Магистратура, 2024

% Схема данных численного эксперимента
% x(1), L(2), x/L(3), dx(4), dx/L(5), x_CO2(6), n_CO2(7), n_AR(8), rho_CO2(9), 
% rho_AR(10), p(11), v(12), T(13), T12(14), T3(15), k(16), a(17), M(18), 
% dVCO2(19), dVAr(20), dCO2(21), dAr(22), dCO2Ar(23), tDCO2(24), tDAr(25),
% lTr(26), lRot(27), lT12(28), lT3(29), trQ(30), vQ12(31), vQ3(32), dQ(33), 
% tDQ(34), sVisc(35), bVisc(36), xxP(37), eFull(38), eT12(39), eT3(40)

%% Инициализация

ps = PlotSettings;
ps.SetFont();

%%
ps.PlotFigs({{'HLLC_bs-100'}, {'Ar'}, {'5'}, {'100'}, ...
    {'3'}, {'0', '1'}, {'0'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta \ne 0', '\zeta = 100\cdot\eta'}, @ps.PlotFigT);
ps.PlotFigs({{'HLLC_bs-100'}, {'Ar'}, {'5'}, {'100'}, ...
    {'3'}, {'0', '1'}, {'0'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta \ne 0', '\zeta = 100\cdot\eta'}, @ps.PlotFigP);
ps.PlotFigs({{'HLLC_bs-100'}, {'Ar'}, {'5'}, {'100'}, ...
    {'3'}, {'0', '1'}, {'0'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta \ne 0', '\zeta = 100\cdot\eta'}, @ps.PlotFigPxx);
ps.PlotFigs({{'HLLC_bs-100'}, {'Ar'}, {'5'}, {'100'}, ...
    {'3'}, {'0', '1'}, {'0'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta \ne 0', '\zeta = 100\cdot\eta'}, @ps.PlotFigQ);
ps.PlotFigs({{'HLLC_bs-100'}, {'Ar'}, {'5'}, {'100'}, ...
    {'3'}, {'0', '1'}, {'0'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta \ne 0', '\zeta = 100\cdot\eta'}, @ps.PlotFigQV);
%%
ps.PlotFigs({{'HLLC_bs-100'}, {'Ar'}, {'5'}, {'100'}, ...
    {'3'}, {'0', '1'}, {'0'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta \ne 0', '\zeta = 100\cdot\eta'}, @ps.PlotFigZE);

%%

ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300', '900'}}, ...
    [10, 11], {'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=900\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=900\,K'}, @ps.PlotFigT);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300'}}, ...
    [10], {'T^{(0)}_{12}=300\,K', 'T^{(0)}_{12}=900\,K'}, @ps.PlotFigT);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigT);
ps.PlotFigs({{'HLLC'}, {'Ar', 'Ne', 'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [2], {'CO_2-Ar', 'CO_2-Ne', 'CO_2-He'}, @ps.PlotFigT);

ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigX);
ps.PlotFigs({{'HLLC'}, {'Ar', 'Ne', 'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [2], {'CO_2-Ar', 'CO_2-Ne', 'CO_2-He'}, @ps.PlotFigX);

ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300'}}, ...
    [10], {'T^{(0)}_{12}=300\,K', 'T^{(0)}_{12}=900\,K'}, @ps.PlotFigP);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigP);
ps.PlotFigs({{'HLLC'}, {'Ar', 'Ne', 'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [2], {'CO_2-Ar', 'CO_2-Ne', 'CO_2-He'}, @ps.PlotFigP);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300'}}, ...
    [10], {'T^{(0)}_{12}=300\,K', 'T^{(0)}_{12}=900\,K'}, @ps.PlotFigPxx);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigPxx);
ps.PlotFigs({{'HLLC'}, {'Ar', 'Ne', 'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [2], {'CO_2-Ar', 'CO_2-Ne', 'CO_2-He'}, @ps.PlotFigPxx);

ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'0', '1'}, {'0', '1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6, 7], {'\zeta=0,\,V_{d}=0', '\zeta=0,\,V_{d}\ne0', ...
    '\zeta\ne0,\,V_{d}=0', '\zeta\ne0,\,V_{d}\ne0'}, @ps.PlotFigPxx);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'50'}, ...
    {'3'}, {'0', '1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta=0', '\zeta\ne0'}, @ps.PlotFigPxx);

ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300'}}, ...
    [10], {'T^{(0)}_{12}=300\,K', 'T^{(0)}_{12}=900\,K'}, @ps.PlotFigQ);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'0', '1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [7], {'V_{d}=0', 'V_{d}\ne0'}, @ps.PlotFigQ);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigQ);
ps.PlotFigs({{'HLLC'}, {'Ar', 'Ne', 'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [2], {'CO_2-Ar', 'CO_2-Ne', 'CO_2-He'}, @ps.PlotFigQ);

ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300'}}, ...
    [10], {'T^{(0)}_{12}=300\,K', 'T^{(0)}_{12}=900\,K'}, @ps.PlotFigQDV);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'0', '1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [7], {'V_{d}=0', 'V_{d}\ne0'}, @ps.PlotFigQDV);
ps.PlotFigs({{'HLLC'}, {'He'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigQDV);
ps.PlotFigs({{'HLLC'}, {'Ar', 'Ne', 'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [2], {'CO_2-Ar', 'CO_2-Ne', 'CO_2-He'}, @ps.PlotFigQDV);

ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300'}}, ...
    [10], {'T^{(0)}_{12}=300\,K', 'T^{(0)}_{12}=900\,K'}, @ps.PlotFigQ);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'50'}, ...
    {'3'}, {'0', '1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [6], {'\zeta=0', '\zeta\ne0'}, @ps.PlotFigQ);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigQ);

ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300', '900'}}, ...
    [10, 11], {'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=900\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=900\,K'}, @ps.PlotFigT);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300'}}, ...
    [10], {'T^{(0)}_{12}=300\,K', 'T^{(0)}_{12}=900\,K'}, @ps.PlotFigP);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300', '900'}}, ...
    [10, 11], {'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=900\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=900\,K'}, @ps.PlotFigQDV);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300', '900'}, {'300', '900'}}, ...
    [10, 11], {'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=300\,K,\,T^{(0)}_{3}=900\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=300\,K', ...
               'T^{(0)}_{12}=900\,K,\,T^{(0)}_{3}=900\,K'}, @ps.PlotFigQV);

ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigT);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigQDV);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigX);

%%
ps.PlotFigs({{'HLLC'}, {'Ar', 'Ne', 'He'}, {'5'}, {'50'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [2], {'CO_2-Ar', 'CO_2-Ne', 'CO_2-He'}, @ps.PlotFigZE);
ps.PlotFigs({{'HLLC'}, {'Ar'}, {'5'}, {'25', '50', '75'}, ...
    {'3'}, {'1'}, {'1'}, {'100'}, {'300'}, {'300'}, {'300'}}, ...
    [4], {'x_{_{CO_2}}=0.25', 'x_{_{CO_2}}=0.50', 'x_{_{CO_2}}=0.75'}, @ps.PlotFigZE);
