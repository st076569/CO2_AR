#ifndef MIXTURECO2AR_H
#define MIXTURECO2AR_H

#include "transport.h"

/// BorderCondition - инструмент расчета макропараметров за ударной волной.
/// Требует информации о характеристиках потока перед УВ, использует для расчета
/// законы сохранения.
class BorderCondition
{
public:

    // Находит равновесные макропараметры за УВ
    void compute(const MacroParam& lP, const int& flag = 1);

    // Возвращает равновесные макропараметры за УВ
    const MacroParam& rP() const;

private:

    // Макропараметры перед УВ и за УВ соответственно
    MacroParam lP_, rP_;

    // Иструмент расчета энергий
    EnergyDc energy_;

    // Общая энергия, массовые доли CO2 и Ar перед УВ соответственно
    double lE_, y0_, y1_, lRho_;

    // Вспомогательные величины
    double alpha_;

    // Уплотнение (1) или разряжение (-1)
    int flag_;

private:

    // Расчет вспомогательных величин
    void initialize(const MacroParam& lP, const int& flag);

    // Возвращает дискриминант
    double computeD(const double& t);

    // Возвращает плотность смеси по температуре
    double computeRho(const double& t);

    // Функция для нахождения равновесной температуры за УВ: F(T_eq) = 0
    double computeF(const double& t);

    // Находит равновесную температуру за УВ методом бисекции
    void computeT();

    // Находит равновесное давление за УВ
    void computeP();

    // Находит равновесные значения плотностей сортов за УВ
    void computeRho();

    // Находит равновесную скорость потока за УВ
    void computeV();
};

/// MixtureCo2Ar - предоставляет инструменты для решения задачи о
/// релаксационных процессах за ударной волной в смеси CO2-Ar
class MixtureCo2Ar
{
public:

    // Шаг по времени, время моделирования
    double dt, time;

    // Хранит номер последней итерации
    int currIter, shockPos;

public:

    MixtureCo2Ar();

    void initialize(const MacroParam& lP, const int& numT = 3,
                    const int& useBVisc = 1, const int& useDiffV = 1);
    void initialize(const QString& name, const double& init_dx);
    void solve(const QString& path, const QString& name, const int& wt);

    QVector<QVector<double>> getMacroParams();

    // Процедура сохранения массива в файл
    void writeToFile(const QString& path, const QString& name,
                     const QVector<QVector<double>>& table);
    void writeToFile(const QString& path, const QString& name,
                     const QVector<double>& vector);

    // Сохранение данных в текстовый файл
    void saveMacroParams(const QString& path, const QString& name,
                         const int& iter);
    void saveTimeLine(const QString& path, const QString& name);

private:

    // Вспомогательные структуры
    QMutex mutex;
    QVector<int> parN_v, parN1_v;
    TemperatureNDc computeT;
    BorderCondition bc;
    int model, useBVisc, useDiffV;

    // Все макропараметры течения во всех точках
    QVector<MacroParam> points;

    // Скорость звука, показатель адиабаты, длины ячеек...
    QVector<double> a, k, dx, freeLength, timeLine;
    QVector<double> trQ, vQ12, vQ3, diffQ, tDiffQ, sVisc, bVisc, xxP;
    QVector<double> eFull, eT12, eT3;
    QVector<double> diffVCO2, diffVAr, diffCO2, diffAr, diffCO2Ar, tDiffCO2, tDiffAr;
    QVector<double> lambdaRot, lambdaTr, lambdaT12, lambdaT3;

    // Консервативные переменные, поточные члены, релаксационные члены
    QVector<QVector<double>> U, F, cF, dF, R, tempU;

    // Значения потока на границах ячеек по методу HLLE
    QVector<QVector<double>> cHllF, dHllF;

private:

    // Пересчет краевого условия на нетеплопроводной стенке
    void heatCondWall(const int& wt);

    // Переходим в систему координат, связанную с трубой
    void transformToRigidCS();

    // Расчет вектора потоков во всех ячейках
    void computeF();

    // Расчет релаксационных членов
    void computeR();

    // Расчет потоков на стыках ячеек методом HLLE
    void computeHlleF();

    // Расчет потоков на стыках ячеек методом HLLС
    void computeHllcF();

    // Расчет потоков на стыках ячеек методом HLL
    void computeHllF();

    // Производит один шаг итерационного процесса
    void step();

    // Возврат к макропараметрам Ui -> points
    void updateMacroParam();

    // Консервативным переменным points -> Ui
    void updateU();

    // Расчет показателей адиабаты 'k' и скоростей звука 'a'
    void updateAK();

    // Рассчитывает временной шаг по критерию Куранта-Фридрихса-Леви
    void updateDt();

    // Находит положение наибольшего разрыва
    void findShockPos();

    // Процедура фильтрации
    void filtrate();

    // Рассчитывает длину свободного пробега в каждой точке
    void updateFreeLength();

    // Обновить правое граничное условие
    void updateRightBC();
};

/// StaticCellCO2Ar - предоставляет инструменты для решения задачи о
/// релаксационных процессах в статической адиабатической или
/// изотермической ячейках
class StaticCellCO2Ar
{
public:

    // Шаг по времени, время моделирования
    double dt, time;

    // Хранит номер последней итерации
    int currIter;

public:

    StaticCellCO2Ar();

    void initialize(const MacroParam& point, const double& dt);
    void initialize(const QString& name);
    void solve();

    QVector<QVector<double>> saveMacroParams();

private:

    // Вспомогательные структуры
    TemperatureNDc computeT;

    // Все макропараметры изотермической и адиабатической ячеек
    MacroParam pointT, pointE;

    // Консервативные переменные
    QVector<double> UT, UE, RT, RE;

    // Переменные для записи в файл
    QVector<double> pT, TT, T12T, T3T;
    QVector<double> pE, TE, T12E, T3E;

private:

    // Расчет релаксационных членов
    void computeR();
    QVector<double> computeR(const MacroParam& point);

    // Производит один шаг итерационного процесса
    void step();

    // Возврат к макропараметрам Ui -> points
    void updateMacroParam();

    // Обновляет значения во временном масштабе
    void updateTimeScale();
};

/// StaticCellCO2ArRelTime - предоставляет инструменты для
/// решения задачи о релаксационных процессах в статической
/// изотермической ячейке
class StaticCellCO2ArRelTime
{
public:

    // Шаг по времени, время моделирования
    double dt, time;

    // Время общей релаксации (сек) * давление (атм)
    double TauP;

public:

    StaticCellCO2ArRelTime();
    void initialize(const MacroParam& point,
                    const double& dt, const QString &name);
    void solve();

    QVector<QVector<double>> saveMacroParams();

private:

    // Вспомогательные структуры
    TemperatureNDc computeT;

    // Все макропараметры изотермической и адиабатической ячеек
    MacroParam pointT;

    // Консервативные переменные
    QVector<double> UT, RT;

    // Текущая полная колебательная энергия
    double vE, vEMax;

private:

    // Расчет релаксационных членов
    void computeR();
    QVector<double> computeR(const MacroParam& point);

    // Производит один шаг итерационного процесса
    void step();

    // Возврат к макропараметрам Ui -> points
    void updateMacroParam();
};

#endif // MIXTURECO2AR_H
