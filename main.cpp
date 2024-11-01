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

    double T = 0.01; // шаг дискретизации (10 мс)
    double g = 9.81; // ускорение свободного падения
    double l = 1.0;  // длина струны (1 м)
    double q = 0.001; // дисперсия шума системы

    double theta_true = 0.1; // начальный угол (радианы)
    double omega_true = 0.0; // начальная угловая скорость

    Visualize visualizer("Pendulum Simulation");
    visualizer.init();

    std::vector<my_Vector> trace;
    visualizer.addTrace(trace);

    int steps = 1000;
    for (int k = 0; k < steps; ++k) {
        double time = k * T;

        my_Vector state(3);
        state.set(0, time);
        state.set(1, theta_true);
        state.set(2, omega_true);

        visualizer.extendsTraceByVec(state);

        // double theta_prev = theta_true;
        // double omega_prev = omega_true;
        //
        // double theta_dot = omega_prev;
        // double omega_dot = - (g / l) * sin(theta_prev);
        //
        // theta_true = theta_prev + theta_dot * T;
        // omega_true = omega_prev + omega_dot * T;

        double theta_n = theta_true;
        double omega_n = omega_true;

        double k1_theta = omega_n;
        double k1_omega = - (g / l) * sin(theta_n);

        double theta_k2 = theta_n + 0.5 * k1_theta * T;
        double omega_k2 = omega_n + 0.5 * k1_omega * T;
        double k2_theta = omega_k2;
        double k2_omega = - (g / l) * sin(theta_k2);

        double theta_k3 = theta_n + 0.5 * k2_theta * T;
        double omega_k3 = omega_n + 0.5 * k2_omega * T;
        double k3_theta = omega_k3;
        double k3_omega = - (g / l) * sin(theta_k3);

        double theta_k4 = theta_n + k3_theta * T;
        double omega_k4 = omega_n + k3_omega * T;
        double k4_theta = omega_k4;
        double k4_omega = - (g / l) * sin(theta_k4);

        theta_true = theta_n + (T / 6.0) * (k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
        omega_true = omega_n + (T / 6.0) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega);


        theta_true += sqrt(q) * randn() * T;
        omega_true += sqrt(q) * randn() * T;
    }

    visualizer.saveTraceToFile("./pendulum_data.txt");
    std::vector<std::pair<double, double>> theta_data;
    std::vector<std::pair<double, double>> omega_data;

    const auto& trace_data = visualizer.getTraces()[0];
    for (const auto& vec : trace_data) {
        double time = vec.get(0);
        double theta = vec.get(1);
        double omega = vec.get(2);

        theta_data.emplace_back(time, theta);
        omega_data.emplace_back(time, omega);
    }

    Gnuplot gp;

    // График угла
    gp << "set terminal png size 800,600\n";
    gp << "set output 'theta.png'\n";
    gp << "set title 'Pendulum Angle Comparison'\n";
    gp << "set xlabel 'Time (s)'\n";
    gp << "set ylabel 'Theta (rad)'\n";
    gp << "plot '-' with lines title 'RK4 Method'\n";
    gp.send1d(theta_data);

    // График угловой скорости
    gp << "set output 'omega.png'\n";
    gp << "set title 'Pendulum Angular Velocity Comparison'\n";
    gp << "set xlabel 'Time (s)'\n";
    gp << "set ylabel 'Omega (rad/s)'\n";
    gp << "plot '-' with lines title 'RK4 Method'\n";;
    gp.send1d(omega_data);

    return 0;
}
