#include "MotionIntegrator.h"
#include <cmath>

struct Parameters {
    double omega;
    double A;
    double B;
};

class GaussNewtonEstimator {
    static double computeVelocity(const Parameters &p, double t) {
        return -p.A * p.omega * sin(p.omega * t) + p.B * p.omega * cos(p.omega * t);
    }

    static void computeJacobian(const Parameters &p, double t,
                                double &dv_dw, double &dv_dA, double &dv_dB) {
        dv_dw = -p.A * (sin(p.omega * t) + p.omega * t * cos(p.omega * t)) +
                p.B * (cos(p.omega * t) - p.omega * t * sin(p.omega * t));
        dv_dA = -p.omega * sin(p.omega * t);
        dv_dB = p.omega * cos(p.omega * t);
    }

    static double computeError(const std::vector<double> &times,
                               const std::vector<double> &velocities,
                               const Parameters &p) {
        double error = 0.0;
        for (size_t i = 0; i < times.size(); i++) {
            double diff = velocities[i] - computeVelocity(p, times[i]);
            error += diff * diff;
        }
        return sqrt(error / times.size());
    }

    static double estimateBaseFrequency(const std::vector<double> &times,
                                        const std::vector<double> &velocities) {
        double mean = 0;
        for (auto v: velocities) mean += v;
        mean /= velocities.size();

        std::vector<double> zero_crossings;
        for (size_t i = 1; i < velocities.size(); i++) {
            double v1 = velocities[i - 1] - mean;
            double v2 = velocities[i] - mean;
            if (v1 * v2 <= 0 && v1 != v2) {
                double t = times[i - 1] - v1 * (times[i] - times[i - 1]) / (v2 - v1);
                zero_crossings.push_back(t);
            }
        }

        if (zero_crossings.size() >= 2) {
            double avg_period = 0;
            for (size_t i = 1; i < zero_crossings.size(); i++) {
                avg_period += zero_crossings[i] - zero_crossings[i - 1];
            }
            avg_period = 2 * avg_period / (zero_crossings.size() - 1);
            return 2 * M_PI / avg_period;
        }
        return 2.0;
    }

public:
    Parameters estimateParameters(const std::vector<double> &times,
                                  const std::vector<double> &velocities,
                                  Parameters initial_guess) {
        const int max_iterations = 100;
        const double tolerance = 1e-6;
        Parameters current = initial_guess;

        for (int iter = 0; iter < max_iterations; iter++) {
            Matrix J(times.size(), 3);
            my_Vector residuals(times.size());

            for (size_t i = 0; i < times.size(); i++) {
                double t = times[i];
                double measured_v = velocities[i];
                double computed_v = computeVelocity(current, t);

                double dv_dw, dv_dA, dv_dB;
                computeJacobian(current, t, dv_dw, dv_dA, dv_dB);

                J.set(i, 0, dv_dw);
                J.set(i, 1, dv_dA);
                J.set(i, 2, dv_dB);

                residuals.set(i, measured_v - computed_v);
            }

            Matrix JT = J.transpose();
            Matrix JTJ = JT * J;
            my_Vector JTr = JT * residuals;

            my_Vector dx = JTJ.solve_system(JTr);

            double alpha = 0.5; // коэффициент демпфирования
            current.omega += alpha * dx.get(0);
            current.A += alpha * dx.get(1);
            current.B += alpha * dx.get(2);

            if (dx.len() < tolerance) {
                break;
            }
        }

        return current;
    }

    Parameters estimateWithMultipleInitialGuesses(const std::vector<double> &times,
                                                  const std::vector<double> &velocities) {
        double base_omega = estimateBaseFrequency(times, velocities);

        std::vector<Parameters> initial_guesses(5);
        for (int i = 0; i < 5; i++) {
            initial_guesses[i].omega = base_omega * (0.8 + 0.4 * i / 4.0);
            initial_guesses[i].A = 1.0;
            initial_guesses[i].B = 0.5;
        }

        Parameters best_estimate = initial_guesses[0];
        double best_error = computeError(times, velocities, best_estimate);

        for (const auto &guess: initial_guesses) {
            Parameters estimate = estimateParameters(times, velocities, guess);
            double error = computeError(times, velocities, estimate);

            if (error < best_error) {
                best_error = error;
                best_estimate = estimate;
            }
        }

        return best_estimate;
    }
};

MotionIntegrator::MotionIntegrator(double step_size, double start_time, double end_time,
                                   double omega, double a_coef, double b_coef)
    : dt(step_size), t_start(start_time), t_end(end_time),
      Omega(omega), A(a_coef), B(b_coef) {
}

MotionIntegrator::MotionIntegrator(double step_size, double start_time, double end_time,
                                   const std::vector<double> &times,
                                   const std::vector<double> &measurements)
    : dt(step_size), t_start(start_time), t_end(end_time) {
    estimateParameters(times, measurements);
}

void MotionIntegrator::estimateParameters(const std::vector<double> &times,
                                          const std::vector<double> &measurements) {
    GaussNewtonEstimator estimator;
    Parameters best_params = estimator.estimateWithMultipleInitialGuesses(times, measurements);
    Omega = best_params.omega;
    A = best_params.A;
    B = best_params.B;
}

my_Vector MotionIntegrator::computeDerivatives(const my_Vector &state, double t) const {
    my_Vector derivatives(2);
    // d(theta)/dt = omega
    derivatives.set(0, state.get(1));
    // d(omega)/dt = -Omega^2 * theta
    derivatives.set(1, -Omega * Omega * state.get(0));
    return derivatives;
}

my_Vector MotionIntegrator::rk4Step(const my_Vector &current_state, double t) const {
    my_Vector k1 = computeDerivatives(current_state, t);

    my_Vector state2 = current_state + k1 * (dt / 2.0);
    my_Vector k2 = computeDerivatives(state2, t + dt / 2.0);

    my_Vector state3 = current_state + k2 * (dt / 2.0);
    my_Vector k3 = computeDerivatives(state3, t + dt / 2.0);

    my_Vector state4 = current_state + k3 * dt;
    my_Vector k4 = computeDerivatives(state4, t + dt);

    return current_state + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

my_Vector MotionIntegrator::theoreticalSolution(double t) const {
    my_Vector state(2);
    state.set(0, A * cos(Omega * t) + B * sin(Omega * t)); //position
    state.set(1, -A * Omega * sin(Omega * t) + B * Omega * cos(Omega * t)); //speed
    return state;
}

std::vector<my_Vector> MotionIntegrator::integrate() {
    std::vector<my_Vector> trajectory;
    my_Vector current_state = theoreticalSolution(t_start);
    trajectory.push_back(current_state);
    double current_time = t_start;
    int steps = static_cast<int>((t_end - t_start) / dt);

    for (int i = 0; i < steps; ++i) {
        current_state = rk4Step(current_state, current_time);
        trajectory.push_back(current_state);
        current_time += dt;
    }
    return trajectory;
}

std::vector<my_Vector> MotionIntegrator::getTheoretical() {
    std::vector<my_Vector> trajectory;
    double t = t_start;
    while (t <= t_end) {
        trajectory.push_back(theoreticalSolution(t));
        t += dt;
    }
    return trajectory;
}

double MotionIntegrator::getMeasurementError(const std::vector<double> &times,
                                             const std::vector<double> &measurements) const {
    double error = 0;
    for (size_t i = 0; i < times.size(); ++i) {
        my_Vector theoretical = theoreticalSolution(times[i]);
        double diff = theoretical.get(0) - measurements[i];
        error += diff * diff;
    }
    return sqrt(error / times.size());
}

my_Vector MotionIntegrator::getFullState(const my_Vector &state, double t) {
    my_Vector full_state(3);
    full_state.set(0, state.get(0)); // theta
    full_state.set(1, state.get(1)); // omega
    full_state.set(2, t);
    return full_state;
}
