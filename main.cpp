#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <iostream>
#include "model/gnuplot-iostream.h"

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

    double T = 0.5;    // шаг дискретизации
    double omega_true = 2.0; // истинная угловая скорость
    double R = 45.5;    // дисперсия шума измерений

    double theta_true = 10.0;     // истинный угол
    double theta_est = 0.0;      // оцененный угол (интегрирование измерений)
    double theta_est_no_noise = 0.0; // оцененный угол без шума (для сравнения)

    std::vector<std::pair<double, double>> theta_true_data;
    std::vector<std::pair<double, double>> theta_est_data;
    std::vector<std::pair<double, double>> theta_est_no_noise_data;
    std::vector<std::pair<double, double>> error_data;

    int steps = 25000;
    for (int k = 0; k < steps; ++k) {
        double time = k * T;

        // истинный угол
        theta_true += omega_true * T;

        // "шумный" угол
        double omega_meas = omega_true + sqrt(R) * randn();

        // интерирование шумного измерения для оценки угла
        theta_est += omega_meas * T;
        // интегрирование истинной углоой скорости без шума (для сравнения)
        theta_est_no_noise += omega_true * T;

        double error = theta_est - theta_true;

        theta_true_data.emplace_back(time, theta_true);
        theta_est_data.emplace_back(time, theta_est);
        theta_est_no_noise_data.emplace_back(time, theta_est_no_noise);
        error_data.emplace_back(time, error);
    }

    Gnuplot gp;

    gp << "set terminal png size 800,600\n";
    gp << "set output 'theta_estimation.png'\n";
    gp << "set title 'True Angle vs Estimated Angle'\n";
    gp << "set xlabel 'Time (ms)'\n";
    gp << "set ylabel 'Theta (rad)'\n";
    gp << "plot '-' with lines title 'True Theta', '-' with lines title 'Estimated Theta', '-' with lines title 'Estimated Theta (No Noise)'\n";
    gp.send1d(theta_true_data);
    gp.send1d(theta_est_data);
    gp.send1d(theta_est_no_noise_data);

    gp << "set output 'error_accumulation.png'\n";
    gp << "set title 'Error Accumulation in Angle Estimation'\n";
    gp << "set xlabel 'Time (ms)'\n";
    gp << "set ylabel 'Error (rad)'\n";
    gp << "plot '-' with lines title 'Estimation Error'\n";
    gp.send1d(error_data);

    return 0;
}
