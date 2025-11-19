#include "ParameterEstimator.h"
#include <cmath>
#include <numeric>

ParameterEstimator::ParameterEstimator(const std::vector<double>& times,
                                     const std::vector<double>& measurements)
    : times(times), velocity_measurements(measurements) {
    double omega_init = estimateInitialFrequency();
    auto [A_init, B_init] = estimateInitialAmplitudes(omega_init);

    parameters = my_Vector(3);
    parameters.set(0, A_init);
    parameters.set(1, B_init);
    parameters.set(2, omega_init);
}

Matrix ParameterEstimator::computeJacobian(const my_Vector& params) const {
    double A = params.get(0);
    double B = params.get(1);
    double Omega = params.get(2);

    Matrix J(times.size(), 3);

    for (size_t i = 0; i < times.size(); ++i) {
        double t = times[i];
        J.set(i, 0, -Omega * std::sin(Omega * t));
        J.set(i, 1, Omega * std::cos(Omega * t));
        J.set(i, 2, -A * (std::sin(Omega * t) + Omega * t * std::cos(Omega * t)) +
                     B * (std::cos(Omega * t) - Omega * t * std::sin(Omega * t)));
    }

    return J;
}

double ParameterEstimator::estimateInitialFrequency() {
    const int max_lag = std::min(100, static_cast<int>(times.size() / 2));
    std::vector<double> autocorr(max_lag);

    double mean = std::accumulate(velocity_measurements.begin(),
                                velocity_measurements.end(), 0.0) / velocity_measurements.size();

    for (int lag = 0; lag < max_lag; lag++) {
        double sum = 0.0;
        for (size_t i = 0; i < velocity_measurements.size() - lag; i++) {
            sum += (velocity_measurements[i] - mean) *
                   (velocity_measurements[i + lag] - mean);
        }
        autocorr[lag] = sum;
    }

    int period = 0;
    for (int i = 1; i < max_lag - 1; i++) {
        if (autocorr[i] > autocorr[i-1] && autocorr[i] > autocorr[i+1]) {
            period = i;
            break;
        }
    }

    return period > 0 ? 2 * M_PI / (period * (times[1] - times[0])) : 1.0;
}

std::pair<double, double> ParameterEstimator::estimateInitialAmplitudes(double omega) {
    Matrix X(times.size(), 2);
    my_Vector y(times.size());

    for (size_t i = 0; i < times.size(); i++) {
        double t = times[i];
        X.set(i, 0, -std::sin(omega * t));
        X.set(i, 1, std::cos(omega * t));
        y.set(i, velocity_measurements[i] / omega);
    }

    Matrix XtX = X.transpose() * X;
    my_Vector Xty = X.transpose() * y;
    my_Vector beta = XtX.solve_system(Xty);

    return {beta.get(0), beta.get(1)};
}

my_Vector ParameterEstimator::computeResiduals(const my_Vector& params) const {
    my_Vector residuals(times.size());

    for (size_t i = 0; i < times.size(); ++i) {
        double A = params.get(0);
        double B = params.get(1);
        double Omega = params.get(2);

        double predicted = -A * Omega * std::sin(Omega * times[i]) +
                          B * Omega * std::cos(Omega * times[i]);
        residuals.set(i, velocity_measurements[i] - predicted);
    }

    return residuals;
}

my_Vector ParameterEstimator::estimate(double tolerance, int max_iterations) {
    my_Vector params = parameters;
    double prev_error = 1e10;

    for (int iter = 0; iter < max_iterations; ++iter) {
        Matrix J = computeJacobian(params);
        my_Vector r = computeResiduals(params);

        double current_error = r.len();
        if (current_error < tolerance) break;

        Matrix JtJ = J.transpose() * J;
        my_Vector Jtr = J.transpose() * r;
        my_Vector delta = JtJ.solve_system(-1.0 * Jtr);

        double alpha = 1.0;
        while (alpha > 1e-10) {
            my_Vector new_params = params + delta * alpha;
            my_Vector new_r = computeResiduals(new_params);
            if (new_r.len() < current_error) {
                params = new_params;
                break;
            }
            alpha *= 0.5;
        }

        if (std::abs(current_error - prev_error) < tolerance)
            break;
        prev_error = current_error;
    }

    return params;
}