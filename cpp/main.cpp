#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <fstream>
#include "model/Visualize.h"

std::vector<double> generateVelocityMeasurements(const std::vector<double>& times,
                                                double omega_true,
                                                double A_true,
                                                double B_true,
                                                double noise_std) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution noise(0.0, noise_std);

    std::vector<double> velocities;
    velocities.reserve(times.size());

    for (double t : times) {
        double true_velocity = -A_true * omega_true * sin(omega_true * t) +
                              B_true * omega_true * cos(omega_true * t);
        velocities.push_back(true_velocity + noise(gen));
    }
    return velocities;
}

void printParameters(const std::string& label, double omega, double A, double B) {
    std::cout << label << ":\n"
              << std::fixed << std::setprecision(6)
              << "Omega = " << omega << "\n"
              << "A     = " << A << "\n"
              << "B     = " << B << "\n" << std::endl;
}

int main() {
    const double t_start = 0.0;
    const double t_end = 100.0;
    const double dt = 0.1;

    const double omega_true = 4.0;
    const double A_true = 2.0;
    const double B_true = 2.0;
    const double noise_std = 0.5;

    std::vector<double> times;
    times.reserve(static_cast<size_t>((t_end - t_start) / dt) + 1);
    for (double t = t_start; t <= t_end; t += dt) {
        times.push_back(t);
    }

    std::vector<double> measurements = generateVelocityMeasurements(
        times, omega_true, A_true, B_true, noise_std
    );


    printParameters("Истинные параметры", omega_true, A_true, B_true);

    std::ofstream iterations_file("integration_results.txt");
    if (iterations_file.is_open()) {
        iterations_file << "Time\tTheta_num\tOmega_num\tTheta_theo\tOmega_theo\tMeasured_v\n";

        iterations_file.close();
    }

    Visualize visualizer("motion_comparison");
    visualizer.init();

    std::vector<my_Vector> numerical_vis;
    std::vector<my_Vector> theoretical_vis;
    std::vector<my_Vector> measurements_vis;


    for (size_t i = 0; i < measurements.size(); ++i) {
        my_Vector measurement_state(3);
        measurement_state.set(0, measurements[i]); // скорость
        measurement_state.set(1, 0.0);            // омега
        measurement_state.set(2, times[i]);       // время
        measurements_vis.push_back(measurement_state);
    }

    visualizer.addTrace(numerical_vis);
    visualizer.addTrace(theoretical_vis);
    visualizer.addTrace(measurements_vis);
    visualizer.saveTraceToFile("motion_results.txt");


    return 0;
}