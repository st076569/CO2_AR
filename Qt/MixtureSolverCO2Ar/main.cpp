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
             const int& nT, const double& p, const double& T)
{
    // Предварительная инициализация солвера
    MixtureCo2Ar solver;
    QString temp, str;
    solver.initialize("energy_data.dat");

    // Вычисление v = M * a
    MacroParam mp(p, 0.0, T, x_CO2);
    EnergyDc e;
    e.initialize();
    e.compute(mp);
    double rho = mp.rho[0] + mp.rho[1];
    double k = 1.0 + mp.p / (rho * e.fullE());
    mp.v = m * qSqrt(k * mp.p / rho);

    // Моделирование течения
    solver.initialize(mp, nT, bVisc, diffV);
    solver.solve();

    // Запись данных в файл
    temp = name;
    temp.append("_M");
    temp.append(str.setNum(qRound(m)));
    temp.append("_x");
    temp.append(str.setNum(qRound(x_CO2 * 100)));
    temp.append("_nT");
    temp.append(str.setNum(nT));
    temp.append("_bV");
    temp.append(str.setNum(bVisc));
    temp.append("_dV");
    temp.append(str.setNum(diffV));
    temp.append("_p");
    temp.append(str.setNum(qRound(p)));
    temp.append("_T");
    temp.append(str.setNum(qRound(T)));
    temp.append(".txt");
    writeToFile(path, temp, solver.saveMacroParams());

    // Вывод сигнальной информации
    std::cout << "\n\n > File name   [-]  : " << temp.toStdString();
    std::cout << "\n > dt          [ms] : " << solver.dt * 1e3;
    std::cout << "\n > Iterations  [-]  : " << solver.currIter << " / "
              << MAX_ITERATION_N << "\n\n";
}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // Путь сохранения результатов
    QString path = "./results";

    // path, name, M, x_CO2, diffV, bVisc, numT, p, T
    //compute(path, "Z_AR", 3.0, 1e-6, 0, 0, 3, 6.66, 300);
    compute(path, "CO2", 5.0, 0.999999, 0, 1, 3, 6.66, 300);

    return 0;
}
