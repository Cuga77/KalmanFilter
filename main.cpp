#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <iostream>
#include "model/Visualize.h"
#include "model/vector.h"
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

    double T = 0.01;    // шаг дискретизации (10 мс)
    double g = 9.81;    // ускорение свободного падения
    double l = 1.0;     // длина маятника (1 м)
    double q = 0.001;   // дисперсия шума системы
    double c = 0.1;     // коэффициент демпфирования :)

    double theta_true = 0.1; //начальный угол (радианы)
    double omega_true = 0.0; //начальная угловая скорость

    std::vector<std::pair<double, double>> phase_data;

    int steps = 15000;
    for (int k = 0; k < steps; ++k) {
        double time = k * T;
        phase_data.emplace_back(theta_true, omega_true);

        double theta_n = theta_true;
        double omega_n = omega_true;

        double k1_theta = omega_n;
        double k1_omega = - (g / l) * sin(theta_n) - c * omega_n;

        double theta_k2 = theta_n + 0.5 * k1_theta * T;
        double omega_k2 = omega_n + 0.5 * k1_omega * T;
        double k2_theta = omega_k2;
        double k2_omega = - (g / l) * sin(theta_k2) - c * omega_k2;

        double theta_k3 = theta_n + 0.5 * k2_theta * T;
        double omega_k3 = omega_n + 0.5 * k2_omega * T;
        double k3_theta = omega_k3;
        double k3_omega = - (g / l) * sin(theta_k3) - c * omega_k3;

        double theta_k4 = theta_n + k3_theta * T;
        double omega_k4 = omega_n + k3_omega * T;
        double k4_theta = omega_k4;
        double k4_omega = - (g / l) * sin(theta_k4) - c * omega_k4;

        theta_true = theta_n + (T / 6.0) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
        omega_true = omega_n + (T / 6.0) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega);

        theta_true += sqrt(q) * randn() * T;
        omega_true += sqrt(q) * randn() * T;
    }

    Gnuplot gp;

    // фазовый портрет
    gp << "set terminal png size 800,600\n";
    gp << "set output 'phase_portrait.png'\n";
    gp << "set title 'Phase Portrait of Damped Pendulum'\n";
    gp << "set xlabel 'Theta (rad)'\n";
    gp << "set ylabel 'Omega (rad/s)'\n";
    gp << "plot '-' with lines title 'Phase Trajectory'\n";
    gp.send1d(phase_data);

    // график угла со временем
    gp << "set output 'theta_damped.png'\n";
    gp << "set title 'Damped Pendulum Angle over Time'\n";
    gp << "set xlabel 'Time (s)'\n";
    gp << "set ylabel 'Theta (rad)'\n";
    gp << "plot '-' with lines title 'Theta (rad)'\n";

    // подготовка данных для графика угла
    std::vector<std::pair<double, double>> theta_time_data;
    for (int k = 0; k < steps; ++k) {
        double time = k * T;
        theta_time_data.emplace_back(time, phase_data[k].first);
    }
    gp.send1d(theta_time_data);

    // график угловой скорости со временем
    gp << "set output 'omega_damped.png'\n";
    gp << "set title 'Damped Pendulum Angular Velocity over Time'\n";
    gp << "set xlabel 'Time (s)'\n";
    gp << "set ylabel 'Omega (rad/s)'\n";
    gp << "plot '-' with lines title 'Omega (rad/s)'\n";

    // подготовка данных для графика угловой скорости
    std::vector<std::pair<double, double>> omega_time_data;
    for (int k = 0; k < steps; ++k) {
        double time = k * T;
        omega_time_data.emplace_back(time, phase_data[k].second);
    }
    gp.send1d(omega_time_data);

    return 0;
}
