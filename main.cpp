#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "model/MotionIntegrator.h"
#include "model/Visualize.h"

std::vector<double> generateMeasurements(const std::vector<double> &times,
                                         double omega_true,
                                         double A_true,
                                         double B_true,
                                         double noise_std) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution noise(0.0, noise_std);
    std::vector<double> measurements;
    for (double t: times) {
        double true_value = A_true * cos(omega_true * t) + B_true * sin(omega_true * t);
        measurements.push_back(true_value + noise(gen));
    }
    return measurements;
}

int main() {
    double t_start = 0.0;
    double t_end = 10.0;
    double dt = 0.1;

    double omega_true = 2.0;
    double A_true = 1.0;
    double B_true = 0.5;
    double noise_std = 0.1;

    std::vector<double> times;
    for (double t = t_start; t <= t_end; t += dt) {
        times.push_back(t);
    }

    std::vector<double> measurements = generateMeasurements(times, omega_true, A_true, B_true, noise_std);
    MotionIntegrator integrator(dt, t_start, t_end, times, measurements);

    double omega_est = integrator.getOmega();
    double A_est = integrator.getACoef();
    double B_est = integrator.getBCoef();

    std::cout << "Истинные параметры:" << std::endl;
    std::cout << "Omega = " << omega_true << std::endl;
    std::cout << "A = " << A_true << std::endl;
    std::cout << "B = " << B_true << std::endl;

    std::cout << "\nОцененные параметры:" << std::endl;
    std::cout << "Omega = " << omega_est << std::endl;
    std::cout << "A = " << A_est << std::endl;
    std::cout << "B = " << B_est << std::endl;

    std::vector<my_Vector> numerical = integrator.integrate();
    std::vector<my_Vector> theoretical = integrator.getTheoretical();

    std::cout << '\n'<< std::endl;
    Visualize visualizer("motion_comparison");
    visualizer.init();

    std::vector<my_Vector> numerical_vis;
    std::vector<my_Vector> theoretical_vis;
    std::vector<my_Vector> measurements_vis;

    double t = t_start;
    for (size_t i = 0; i < numerical.size(); ++i) {
        numerical_vis.push_back(MotionIntegrator::getFullState(numerical[i], t));
        theoretical_vis.push_back(MotionIntegrator::getFullState(theoretical[i], t));
        t += dt;
    }

    for (size_t i = 0; i < measurements.size(); ++i) {
        my_Vector measurement_state(3);
        measurement_state.set(0, measurements[i]); // theta
        measurement_state.set(1, 0.0); // omega
        measurement_state.set(2, times[i]);
        measurements_vis.push_back(measurement_state);
    }

    visualizer.addTrace(numerical_vis);
    visualizer.addTrace(theoretical_vis);
    visualizer.addTrace(measurements_vis);

    visualizer.saveTraceToFile("motion_results.txt");

    double error = integrator.getMeasurementError(times, measurements);
    std::cout << "\nСреднеквадратичная ошибка: " << error << std::endl;

    return 0;
}
