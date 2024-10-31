#ifndef GLOBAL_H
#define GLOBAL_H

#include <QVector>
#include <QtMath>
#include <QtNumeric>
#include <QDebug>
#include <QObject>
#include <QThread>
#include <QMutexLocker>
#include <QFutureWatcher>
#include <QtConcurrent>
#include <iostream>

// Настройка температурного диапазона
#define T_MAX 5000.0
#define T_MIN 200.0
#define T_NUM 9601
#define Y_NUM 401

// Малые величины для оценки точности
#define EPSILON 1e-10

// Число итераций, положение разрыва, число Куранта
#define RELATIVE_SHOCK_POSITION 0.001
#define MAX_ITERATION_N 20000
#define CFL 0.9

// Параметр химической добавки в смеси {0-He', 1-Ne', 2-Ar', 3-Kr, 4-Xe}
#define ADMIXTURE 2

// Параметры сетки решателя
#define GRID_N 101
#define GRID_L 0.5
#define GRID_DELTA_0 2.0
#define GRID_DELTA_1 1.0
#define GRID_DEGREE 5

// Постоянная Больцмана
#define K_BOLTZMANN 1.3805e-23

// Атомная единица массы
#define ATOMIC_MASS_UNIT 1.6605402e-27

// Постоянная Планка и скорость света в вакууме
#define H_PLANCK 6.6254e-34
#define C_LIGHT 2.99792458e8

// Число Авогадро
#define N_AVOGADRO 6.0221e+23

// Размерный коэффициент для расчета энергии
#define W_W 1.60219e-19 / 8065.47

// Кол-во сортов смеси
#define N_SORTS 2
#define N_POLYATOMIC_SORTS 1

// Характеристики внутренних степеней свободы CO2
#define N_VIBR_DEGREES 3
#define N_VIBR_L1 30
#define N_VIBR_L2 65
#define N_VIBR_L3 20

// Энергия диссоциации : CO2
#define DISS_ENERGY 64017.0

// Настройка метода Гаусса
#define N_MAX 30

// Число уравнений системы
#define SYSTEM_ORDER 6

/// Структура данных смеси СO2-Ar:
/// 1) Спектроскопические данные для молекулы СО2;
/// 2) Массы молекул сортов СO2, Ar;
/// 3) Газокинетические диаметры молекул сортов СO2, Ar;
/// 4) Глубина потенциальной ямы для молекул сортов : СO2, Ar.
class Mixture
{
public:

    // Методы расчета вспомогательных величин
    static double mass(const int& i);
    static double sigma(const int& i);
    static double epsilon(const int& i);
    static double vEnergy(const int& i, const int& j, const int& k);
    static double reducedMass(const int& i, const int& j);

    // Времена колебательной релаксации
    static double tauVTCO2CO2(const double& t, const double& rho_CO2);
    static double tauVVCO2CO2(const double& t, const double& rho_CO2);
    static double tauVTCO2Ar(const double& t, const double& rho_Ar);
    static double tauVVCO2Ar(const double& t, const double& rho_Ar);
};

// Основные макропараметры (rho[0] -> CO2, rho[1] -> Ar)
class MacroParam
{
public:

    // Основные макропараметры (не независимые)
    QVector<double> rho;
    double p, v, t, t12, t3;

public:


    // Инициализация по умолчанию и для равновесного и неравновесного случаев
    MacroParam();
    MacroParam(const double& p, const double& v, const double& t,
               const double& x_CO2);
    MacroParam(const double& p, const double& v, const double& t,
               const double& t12, const double& t3, const double& x_CO2);

    // Инициализация плотностей, зная {p, t, x_CO2 in [0, 1]}
    void computeRho(const double& x_CO2);

    // Инициализация для равновесного случая
    void initialize(const double& p, const double& v, const double& t,
                    const double& x_CO2);
};

/// ProgressBar - реализует шкалу прогресса в консоли
class ProgressBar
{
public:

    // Текущее и предыдущее заполнение, шаг
    double n0, n1, dx;

public:

    // Обнуление
    ProgressBar();

    // Инициализация для равновесного случая
    void initialize(const double& t);

    // Обновить шкалу
    void update(const double& dt);
};

#endif // GLOBAL_H
