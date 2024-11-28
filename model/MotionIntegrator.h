#pragma once
#include <vector>
#include <fstream>
#include "vector.h"
#include "matrix.h"

struct Parameters {
    double omega;
    double A;
    double B;

    Parameters() : omega(0), A(0), B(0) {}
    Parameters(double w, double a, double b) : omega(w), A(a), B(b) {}
};

class MotionIntegrator {
private:
    double dt;
    double t_start;
    double t_end;

    double Omega;
    double A;
    double B;

    std::ofstream log_file;

    my_Vector computeDerivatives(const my_Vector& state, double t) const;
    my_Vector rk4Step(const my_Vector& current_state, double t) const;
    my_Vector theoreticalSolution(double t) const;

    void logIteration(int iter, double omega, double A, double B, double error);

public:
    MotionIntegrator(double step_size, double start_time, double end_time,
                     double omega, double a_coef, double b_coef);

    MotionIntegrator(double step_size, double start_time, double end_time,
                     const std::vector<double>& times,
                     const std::vector<double>& measurements);

    ~MotionIntegrator() {
        if(log_file.is_open()) {
            log_file.close();
        }
    }

    void estimateParameters(const std::vector<double>& times,
                          const std::vector<double>& measurements);

    std::vector<my_Vector> integrate();
    std::vector<my_Vector> getTheoretical();

    double getMeasurementError(const std::vector<double>& times,
                             const std::vector<double>& measurements) const;

    double getOmega() const { return Omega; }
    double getACoef() const { return A; }
    double getBCoef() const { return B; }

    static my_Vector getFullState(const my_Vector& state, double t);
};