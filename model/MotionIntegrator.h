#pragma once
#include <vector>
#include "vector.h"
#include "matrix.h"

class MotionIntegrator {
    double dt;
    double t_start;
    double t_end;

    double Omega;
    double A; // амплитуда при косинусе
    double B; // амплитуда при синусе

    my_Vector computeDerivatives(const my_Vector &state, double t) const;

    my_Vector rk4Step(const my_Vector &current_state, double t) const;

    my_Vector theoreticalSolution(double t) const;

    double estimateInitialOmega(const std::vector<double> &times,
                                const std::vector<double> &measurements);

public:
    MotionIntegrator(double step_size, double start_time, double end_time,
                     double omega, double a_coef, double b_coef);

    MotionIntegrator(double step_size, double start_time, double end_time,
                     const std::vector<double> &times,
                     const std::vector<double> &measurements);

    std::vector<my_Vector> integrate();

    std::vector<my_Vector> getTheoretical();

    void estimateParameters(const std::vector<double> &times,
                            const std::vector<double> &measurements);

    double getOmega() const { return Omega; }
    double getACoef() const { return A; }
    double getBCoef() const { return B; }

    static my_Vector getFullState(const my_Vector &state, double t);

    double getMeasurementError(const std::vector<double> &times,
                               const std::vector<double> &measurements) const;
};

struct OptimizationParams {
    const std::vector<double> &times;
    const std::vector<double> &measurements;

    OptimizationParams(const std::vector<double> &t, const std::vector<double> &m)
        : times(t), measurements(m) {
    }
};

double errorFunction(const std::vector<double> &params, void *opt_data);
