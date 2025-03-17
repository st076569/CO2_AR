#include <QCoreApplication>
#include <QFile>
#include <QTextStream>
#include <QString>
#include <QRandomGenerator>
#include "mixture.h"

// Запись данных в файл
void writeToFile(const QString& path, const QString& name,
                 const QVector<QVector<double>>& table)
{
    // Открываем файл для записи
    QFile file1(path + "/" + name);
    file1.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out1(&file1);

    // Настройка параметров вывода
    out1.setFieldWidth(16);
    out1.setRealNumberPrecision(8);
    out1.setRealNumberNotation(QTextStream::ScientificNotation);

    // Запись таблицы в файл
    for (int i = 0; i < table[0].size(); ++i)
    {
        for (int j = 0; j < table.size(); ++j)
        {
            out1 << table[j][i];
        }
        out1 << '\n';
    }

    // Закрытие файла
    file1.close();
}

// Процедура реализации варьирования параметров
void compute(const QString& path, const QString& name, const double& m,
             const double& x_CO2, const int& diffV, const int& bVisc,
             const int& nT, const double& p, const double& T,
             const double& T12, const double& T3, const int& wt)
{
    // Предварительная инициализация солвера
    MixtureCo2Ar solver;
    QString temp1, str;

    // Вычисление v = M * a
    MacroParam mp(p, 0.0, T, T12, T3, x_CO2);
    EnergyDc e;
    e.initialize();
    e.compute(mp);

    // Пересчет макропараметров
    double rho = mp.rho[0] + mp.rho[1];
    double k = 1.0 + mp.p / (rho * e.fullE());
    double n0 = mp.rho[0] / Mixture::mass(0);
    double n1 = mp.rho[1] / Mixture::mass(1);
    double x0 = n0 / (n0 + n1);
    double d = x0 * Mixture::sigma(0) + (1.0 - x0) * Mixture::sigma(1);
    double gd = (GRID_DELTA_1 - GRID_DELTA_0) * x0 + GRID_DELTA_0;
    double dx_0 = gd * K_BOLTZMANN * mp.t / qSqrt(2) / M_PI / qPow(d, 2.0) / mp.p;
    mp.v = m * qSqrt(k * mp.p / rho);

    // Формируем имя файла
    temp1 = name;
    switch (ADMIXTURE)
    {
        case 0:
            temp1.append("_He");
        break;
        case 1:
            temp1.append("_Ne");
        break;
        case 2:
            temp1.append("_Ar");
        break;
        case 3:
            temp1.append("_Kr");
        break;
        case 4:
            temp1.append("_Xe");
        break;
        default:
            temp1.append("_Ar");
    }
    temp1.append("_M-");
    temp1.append(str.setNum(qRound(m)));
    temp1.append("_x-");
    temp1.append(str.setNum(qRound(x_CO2 * 100)));
    temp1.append("_nT-");
    temp1.append(str.setNum(nT));
    temp1.append("_bV-");
    temp1.append(str.setNum(bVisc));
    temp1.append("_dV-");
    temp1.append(str.setNum(diffV));
    temp1.append("_p-");
    temp1.append(str.setNum(qRound(p)));
    temp1.append("_T-");
    temp1.append(str.setNum(qRound(T)));
    temp1.append("_T12-");
    temp1.append(str.setNum(qRound(T12)));
    temp1.append("_T3-");
    temp1.append(str.setNum(qRound(T3)));
    temp1.append("_WT-");
    temp1.append(str.setNum(wt));

    // Моделирование течения
    solver.initialize("energy_data_1.dat", dx_0);
    solver.initialize(mp, nT, bVisc, diffV);
    solver.solve(path, temp1, wt);

    // Вывод сигнальной информации
    std::cout << "\n\n > File name  [.csd] : " << temp1.toStdString();
    std::cout << "\n > dt           [us] : " << solver.dt * 1e6;
    std::cout << "\n > Time         [us] : " << solver.time * 1e6;
    std::cout << "\n\n";
}

/*
// Процедура реализации варьирования параметров
void computeVar(const QString& path, const QString& name,
                const QVector<double>& m, const QVector<double>& x_CO2,
                const QVector<int>& diffV, const QVector<int>& bVisc,
                const QVector<int>& nT, const QVector<double>& p,
                const QVector<double>& T, const QVector<double>& T12,
                const QVector<double>& T3)
{
    for (int i1 = 0; i1 < m.size();     ++i1)
    for (int i2 = 0; i2 < x_CO2.size(); ++i2)
    for (int i3 = 0; i3 < diffV.size(); ++i3)
    for (int i4 = 0; i4 < bVisc.size(); ++i4)
    for (int i5 = 0; i5 < nT.size();    ++i5)
    for (int i6 = 0; i6 < p.size();     ++i6)
    for (int i7 = 0; i7 < T.size();     ++i7)
    for (int i8 = 0; i8 < T12.size();   ++i8)
    for (int i9 = 0; i9 < T3.size();    ++i9)
        compute(path, name, m[i1], x_CO2[i2], diffV[i3], bVisc[i4], nT[i5],
                  p[i6], T[i7], T12[i8], T3[i9]);
}

// Процедура реализации варьирования параметров
void computeSC(const QString& path, const QString& name, const double& x_CO2,
             const double& p, const double& T, const double& T12,
             const double& T3, const double& dt)
{
    // Предварительная инициализация солвера
    StaticCellCO2Ar solver;
    QString temp, str;
    solver.initialize("energy_data_1.dat");
    MacroParam mp(p, 0.0, T, T12, T3, x_CO2);

    // Моделирование течения
    solver.initialize(mp, dt);
    solver.solve();

    // Запись данных в файл
    temp = name;
    temp.append("_x");
    temp.append(str.setNum(qRound(x_CO2 * 100)));
    temp.append("_p");
    temp.append(str.setNum(qRound(p)));
    temp.append("_T");
    temp.append(str.setNum(qRound(T)));
    temp.append("_12T");
    temp.append(str.setNum(qRound(T12)));
    temp.append("_3T");
    temp.append(str.setNum(qRound(T3)));
    temp.append("_dt");
    temp.append(str.setNum(qRound(dt * 1e9)));
    temp.append(".tsd");
    writeToFile(path, temp, solver.saveMacroParams());

    // Вывод сигнальной информации
    std::cout << "\n\n > File name  [.tsd] : " << temp.toStdString();
    std::cout << "\n > dt           [ns] : " << solver.dt * 1e9;
    std::cout << "\n > Time         [us] : " << solver.time * 1e6;
    std::cout << "\n > Iterations    [-] : " << solver.currIter << " / "
              << MAX_ITERATION_N << "\n\n";
}
*/

// Процедура реализации варьирования параметров
void computeSCRelTime(const QString& path, const QString& name,
                      const double& x_CO2, const double& p,
                      const QVector<double>& T, const double& T12,
                      const double& T3, const double& dt)
{
    QVector<QVector<double>> table;
    table.resize(3);
    QString temp, str;

    for (int i = 0; i < T.size(); ++i)
    {
        StaticCellCO2ArRelTime solver;
        MacroParam mp(p, 0.0, T[i], T12, T3, x_CO2);

        solver.initialize(mp, dt, "energy_data_1.dat");
        solver.solve();

        table[0].append(T[i]);
        table[1].append(qPow(T[i], -1.0 / 3.0));
        table[2].append(solver.TauP);
    }

    // Запись данных в файл
    temp = name;
    switch (ADMIXTURE)
    {
        case 0:
            temp.append("_He");
        break;
        case 1:
            temp.append("_Ne");
        break;
        case 2:
            temp.append("_Ar");
        break;
        case 3:
            temp.append("_Kr");
        break;
        case 4:
            temp.append("_Xe");
        break;
        default:
            temp.append("_Ar");
    }
    temp.append("_x-");
    temp.append(str.setNum(qRound(x_CO2 * 100)));
    temp.append("_p-");
    temp.append(str.setNum(qRound(p / 101325)));
    temp.append("_T12-");
    temp.append(str.setNum(qRound(T12)));
    temp.append("_T3-");
    temp.append(str.setNum(qRound(T3)));
    temp.append("_dt-");
    temp.append(str.setNum(qRound(dt * 1e9)));
    temp.append(".tsd");
    writeToFile(path, temp, table);

    // Вывод сигнальной информации
    std::cout << "\n\n > File name [.tsd] : " << temp.toStdString() << "\n\n";
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // Путь сохранения результатов
    QString path = "./results_answer";

    // path, name, M, x_CO2, diffV, bVisc, numT, p, T, T12, T3
    compute(path, "HLL", 5.0, 0.5, 1, 1, 3, 100, 300, 300, 300, 0);
    //computeNE(path, "HLLC", 5.0, 0.5, 1, 1, 3, 100, 300, 300, 300);
    //computeVar(path, "HLLC", {5.0}, {0.5}, {0, 1}, {0, 1}, {3}, {100}, {300}, {300}, {300});
    //computeVar(path, "HLLC", {5.0}, {1e-7, 0.25, 0.75, 1 - 1e-7}, {1}, {1}, {3}, {100}, {300},
    //{300}, {300});
    //computeVar(path, "HLLC", {5.0}, {0.9999}, {0}, {1}, {3}, {100}, {300}, {300}, {300});
    //computeVar(path, "HLLC", {5.0}, {1e-6, 0.25, 0.75, 1 - 1e-6}, {1}, {1}, {3}, {100}, {300},
    //{300}, {300});

    // path, name, x_CO2, p, T, T12, T3, dt
    //computeSC(path, "SC", 0.999999, 80, 1800, 2000, 2000, 1e-6);

    //QVector<double> t = {300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
    //                     1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100};
    //computeSCRelTime(path, "SCRT", 0.6, 101325, t, 220, 220, 1e-9);

    // Запись таблицы энергий в файл
    //TemperatureNDc computeT;
    //computeT.initialize(200, 5000, 9601, 401);
    //computeT.writeToFile("energy_data_1.dat");

    return 0;
}
