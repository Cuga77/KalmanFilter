#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <fstream>

#include "model/matrix.h"
#include "model/vector.h"
#include "model/KalmanFilter.h"

// Helper to avoid linker errors with the existing library
Matrix createIdentity(int n) {
    Matrix res(n, n, 0.0);
    for (int i = 0; i < n; ++i) res.set(i, i, 1.0);
    return res;
}

int main() {
    // ==========================================
    // 1. SYSTEM PARAMETERS (Differential Equation setup)
    // ==========================================
    const double t_end = 10.0;
    const double dt = 0.01; 
    const double Omega_true = 3.0; // Natural frequency
    const double Initial_Amplitude = 1.0; 
    const double sigma_meas = 0.2; // Measurement noise (Sensor)

    std::cout << "Starting Simulation..." << std::endl;
    std::cout << "Freq (Omega): " << Omega_true << ", dt: " << dt << std::endl;

    // ==========================================
    // 2. DATA GENERATION (Ground Truth)
    // ==========================================
    std::vector<double> times, true_angles, measured_velocities;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> noise(0.0, sigma_meas);

    for (double t = 0; t <= t_end; t += dt) {
        times.push_back(t);
        // Exact solution of phi'' + w^2*phi = 0
        double phi = Initial_Amplitude * std::cos(Omega_true * t);
        double phi_dot = -Initial_Amplitude * Omega_true * std::sin(Omega_true * t);
        
        true_angles.push_back(phi);
        measured_velocities.push_back(phi_dot + noise(gen));
    }

    // ==========================================
    // 3. KALMAN FILTER SETUP
    // ==========================================
    Matrix x0(2, 1); x0.set(0, 0, Initial_Amplitude); x0.set(1, 0, 0.0);
    Matrix P0 = createIdentity(2);
    Matrix H(1, 2); H.set(0, 0, 0.0); H.set(0, 1, 1.0); // Measuring velocity
    Matrix R(1, 1); R.set(0, 0, sigma_meas * sigma_meas);

    // --- MODEL 1: Constant Velocity (Kinematic) ---
    // ODE: x'' = 0
    // Fundamental Matrix: [[1, dt], [0, 1]] (Taylor series, order 1)
    Matrix F_perm(2, 2);
    F_perm.set(0, 0, 1.0); F_perm.set(0, 1, dt);
    F_perm.set(1, 0, 0.0); F_perm.set(1, 1, 1.0);
    
    // High Process Noise (Q) because the model is physically wrong
    Matrix Q_perm = createIdentity(2);
    Q_perm.set(0, 0, 1e-4); Q_perm.set(1, 1, 0.1);  
    
    KalmanFilter kf_perm(x0, P0, F_perm, Q_perm, R, H);

    // --- MODEL 2: Harmonic Oscillator (Dynamic) ---
    // ODE: x'' + w^2*x = 0
    // Fundamental Matrix: Exact trigonometric solution
    Matrix F_osc(2, 2);
    double c = std::cos(Omega_true * dt);
    double s = std::sin(Omega_true * dt);
    F_osc.set(0, 0, c);                 F_osc.set(0, 1, s / Omega_true);
    F_osc.set(1, 0, -Omega_true * s);   F_osc.set(1, 1, c);

    // Low Process Noise (Q) because the model matches reality
    Matrix Q_osc = createIdentity(2);
    Q_osc.set(0, 0, 1e-6); Q_osc.set(1, 1, 1e-6);

    KalmanFilter kf_osc(x0, P0, F_osc, Q_osc, R, H);

    // ==========================================
    // 4. SIMULATION LOOP
    // ==========================================
    std::ofstream outfile("results.txt");
    outfile << "Time TrueAngle MeasVel AnglePerm AngleOsc" << std::endl;

    double error_sq_perm = 0.0;
    double error_sq_osc = 0.0;

    for (size_t i = 0; i < times.size(); ++i) {
        double z = measured_velocities[i];

        // Predict
        kf_perm.predict();
        kf_osc.predict();

        // Update
        kf_perm.update(z);
        kf_osc.update(z);

        // Errors
        double diff_perm = true_angles[i] - kf_perm.getAngle();
        double diff_osc = true_angles[i] - kf_osc.getAngle();
        error_sq_perm += diff_perm * diff_perm;
        error_sq_osc += diff_osc * diff_osc;

        // Save
        outfile << times[i] << " " << true_angles[i] << " " << z << " "
                << kf_perm.getAngle() << " " << kf_osc.getAngle() << std::endl;
    }
    outfile.close();

    // ==========================================
    // 5. REPORT
    // ==========================================
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "RESULTS (Mean Squared Error Comparison):" << std::endl;
    std::cout << "  Model 1 (Constant Velocity): " << error_sq_perm << std::endl;
    std::cout << "  Model 2 (Harmonic Osc.):     " << error_sq_osc << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    if (error_sq_osc < error_sq_perm) {
        std::cout << "SUCCESS: Physical model outperforms kinematic model." << std::endl;
        std::cout << "Improvement factor: " << (error_sq_perm / error_sq_osc) << "x" << std::endl;
    } else {
        std::cout << "FAIL: Check Q/R parameters." << std::endl;
    }

    return 0;
}