#include <cstdlib>
#include <cmath>
#include <ctime>
#include "model/vector.cpp"
#include "model/matrix.cpp"
#include "model/KalmanFilter.cpp"
#include "model/Visualize.cpp"

double randn() {
    static bool hasSpare = false;
    static double rand1, rand2;

    if (hasSpare)
    {
        hasSpare = false;
        return sqrt(rand1) * sin(rand2);
    }

    hasSpare = true;

    rand1 = rand() / ((double) RAND_MAX);
    if (rand1 < 1e-100) rand1 = 1e-100;
    rand1 = -2 * log(rand1);
    rand2 = (rand() / ((double) RAND_MAX)) * M_PI * 2;

    return sqrt(rand1) * cos(rand2);
}

int main() {
    srand(time(NULL));

    double T = 0.01; // щаг дискретизации (10 мс)
    double g = 9.81; // ускорение свободного падения
    double l = 1.0;  // длина струны (1 м)
    double q = 0.001; // дисперсия шума системы
    double r = 0.01;  // дисперсия шума измерений

    //начальное состояние
    my_Vector x0(2);
    x0.set(0, 0.1); // начальный угол (радианы)
    x0.set(1, 0.0); // начальная угловая скорость

    // Начальная ковариационная матрица ошибки оценки
    Matrix P0(2, 2);
    P0.set(0, 0, 0.1);
    P0.set(0, 1, 0.0);
    P0.set(1, 0, 0.0);
    P0.set(1, 1, 0.1);

    KalmanFilter kf(T, g, l, q, r);
    kf.initialize(x0, P0);

    Visualize visualizer("Pendulum Simulation");
    visualizer.init();

    std::vector<my_Vector> trajectory;

    double omega_natural = sqrt(g / l);
    double theta_true = x0.vec[0];
    double omega_true = x0.vec[1];

    int steps = 1000;
    for (int k = 0; k < steps; ++k) {
        double time = k * T;

        // Обновление истинного состояния (без учета шума системы)
        double theta_prev = theta_true;
        theta_true = theta_true + omega_true * T;
        omega_true = omega_true - (g / l) * theta_prev * T;

        // Добавляем шум системы
        theta_true += sqrt(q) * randn() * T;
        omega_true += sqrt(q) * randn() * T;

        // Генерация измерения с шумом
        double measurement_noise = sqrt(r) * randn();
        double measurement = theta_true + measurement_noise;

        my_Vector z(1);
        z.set(0, measurement);

        kf.predict();
        kf.update(z);

        my_Vector x_est = kf.getState();

        trajectory.push_back(x_est);

        // visualizer.extendsTraceByVec(x_est);
    }

    visualizer.addTrace(trajectory);

    // visualizer.printTraces();

    return 0;
}
