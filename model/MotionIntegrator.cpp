#include "MotionIntegrator.h"
#include <cmath>
#include <fstream>

class GaussNewtonEstimator {
    std::ofstream& log_file;

    static double computeVelocity(const Parameters& p, double t) {
        return -p.A * p.omega * sin(p.omega * t) + p.B * p.omega * cos(p.omega * t);
    }

    static Parameters computeInitialGuess(const std::vector<double>& times,
                                        const std::vector<double>& velocities) {
        Parameters guess;

        int zero_crossings = 0;
        for(size_t i = 1; i < velocities.size(); i++) {
            if(velocities[i-1] * velocities[i] <= 0) {
                zero_crossings++;
            }
        }
        
        double T = times.back() - times.front();
        guess.omega = M_PI * zero_crossings / T;
        
        double max_vel = 0;
        for(auto v : velocities) {
            max_vel = std::max(max_vel, std::abs(v));
        }
        
        guess.A = max_vel / (guess.omega * sqrt(2.0));
        guess.B = guess.A;
        
        return guess;
    }

public:
    GaussNewtonEstimator(std::ofstream& log) : log_file(log) {}

    Parameters estimateParameters(const std::vector<double>& times,
                                const std::vector<double>& velocities) {
        Parameters current = computeInitialGuess(times, velocities);

        const int max_iterations = 1000;
        const double tolerance = 1e-10;
        double prev_error = std::numeric_limits<double>::max();

        for(int iter = 0; iter < max_iterations; iter++) {
            Matrix J(times.size(), 3);
            my_Vector residuals(times.size());

            for(size_t i = 0; i < times.size(); i++) {
                double t = times[i];
                double measured_v = velocities[i];
                double computed_v = computeVelocity(current, t);

                double dv_dw = -current.A * (sin(current.omega * t) + current.omega * t * cos(current.omega * t)) +
                               current.B * (cos(current.omega * t) - current.omega * t * sin(current.omega * t));
                double dv_dA = -current.omega * sin(current.omega * t);
                double dv_dB = current.omega * cos(current.omega * t);

                J.set(i, 0, dv_dw);
                J.set(i, 1, dv_dA);
                J.set(i, 2, dv_dB);

                residuals.set(i, measured_v - computed_v);
            }

            Matrix JT = J.transpose();
            Matrix JTJ = JT * J;
            my_Vector JTr = JT * residuals;

            my_Vector dx = JTJ.solve_system(JTr);

            double alpha = 1.0;
            Parameters next = current;
            double current_error;

            do {
                next.omega = current.omega + alpha * dx.get(0);
                next.A = current.A + alpha * dx.get(1);
                next.B = current.B + alpha * dx.get(2);
                
                next.omega = std::max(0.1, std::min(10.0, next.omega));

                current_error = 0;
                for(size_t i = 0; i < times.size(); i++) {
                    double diff = velocities[i] - computeVelocity(next, times[i]);
                    current_error += diff * diff;
                }
                current_error = sqrt(current_error / times.size());

                alpha *= 0.5;
            } while(current_error >= prev_error && alpha > 1e-15);

            log_file << iter << "\t"
                    << next.omega << "\t"
                    << next.A << "\t"
                    << next.B << "\t"
                    << current_error << "\n";

            if(fabs(current_error - prev_error) < tolerance) {
                break;
            }

            current = next;
            prev_error = current_error;
        }

        return current;
    }
};

MotionIntegrator::MotionIntegrator(double step_size, double start_time, double end_time,
                                   double omega, double a_coef, double b_coef)
    : dt(step_size), t_start(start_time), t_end(end_time),
      Omega(omega), A(a_coef), B(b_coef) {
}

MotionIntegrator::MotionIntegrator(double step_size, double start_time, double end_time,
                                   const std::vector<double>& times,
                                   const std::vector<double>& measurements)
    : dt(step_size), t_start(start_time), t_end(end_time) {
    estimateParameters(times, measurements);
}


void MotionIntegrator::estimateParameters(const std::vector<double>& times,
                                        const std::vector<double>& measurements) {
    log_file.open("parameter_estimation.log");
    if (!log_file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл для логирования оценки параметров");
    }
    
    log_file << "# Iteration\tOmega\tA\tB\tError\tAlpha\n";
    
    GaussNewtonEstimator estimator(log_file);
    Parameters estimated = estimator.estimateParameters(times, measurements);
    
    Omega = estimated.omega;
    A = estimated.A;
    B = estimated.B;
    
    log_file << "\n# Final parameters:\n";
    log_file << "# Omega = " << Omega << "\n";
    log_file << "# A = " << A << "\n";
    log_file << "# B = " << B << "\n";
    log_file << "# Final error = " << getMeasurementError(times, measurements) << "\n";
}

void MotionIntegrator::logIteration(int iter, double omega, double A, double B, double error) {
    if (log_file.is_open()) {
        log_file << iter << "\t" << omega << "\t" << A << "\t" << B << "\t" << error << "\n";
    }
}

my_Vector MotionIntegrator::computeDerivatives(const my_Vector& state, double t) const {
    my_Vector derivatives(2);
    derivatives.set(0, state.get(1));
    derivatives.set(1, -Omega * Omega * state.get(0));
    return derivatives;
}

my_Vector MotionIntegrator::rk4Step(const my_Vector& current_state, double t) const {
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
    state.set(0, A * cos(Omega * t) + B * sin(Omega * t));
    state.set(1, -A * Omega * sin(Omega * t) + B * Omega * cos(Omega * t));
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

double MotionIntegrator::getMeasurementError(const std::vector<double>& times,
                                           const std::vector<double>& measurements) const {
    double error = 0.0;
    for (size_t i = 0; i < times.size(); ++i) {
        my_Vector theoretical = theoreticalSolution(times[i]);
        double diff = theoretical.get(0) - measurements[i];
        error += diff * diff;
    }
    return sqrt(error / times.size());
}

my_Vector MotionIntegrator::getFullState(const my_Vector& state, double t) {
    my_Vector full_state(3);
    full_state.set(0, state.get(0));
    full_state.set(1, state.get(1));
    full_state.set(2, t);
    return full_state;
}