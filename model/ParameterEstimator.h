#pragma once

#include "matrix.h"
#include "vector.h"
#include <vector>
#include <utility>

class ParameterEstimator {
    std::vector<double> times;
    std::vector<double> velocity_measurements;
    my_Vector parameters; // [A, B, Omega]

    Matrix computeJacobian(const my_Vector& params) const;
    my_Vector computeResiduals(const my_Vector& params) const;

    double estimateInitialFrequency();
    std::pair<double, double> estimateInitialAmplitudes(double omega);

public:
    ParameterEstimator(const std::vector<double>& times,
                      const std::vector<double>& measurements);
    my_Vector estimate(double tolerance = 1e-6, int max_iterations = 100);
};