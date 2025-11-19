#include "KalmanFilter.h"
#include <cmath>

KalmanFilter::KalmanFilter(Matrix x0, Matrix P0, Matrix F_in, Matrix Q_in, Matrix R_in, Matrix H_in)
    : x(x0), P(P0), F(F_in), Q(Q_in), R(R_in), H(H_in)
{
}

void KalmanFilter::setF(const Matrix& newF) {
    F = newF;
}

void KalmanFilter::predict() {
    // Prediction Step based on Discrete Linear Dynamic System:
    // x(k+1) = F * x(k)
    // P(k+1) = F * P(k) * F^T + Q

    // 1. State prediction (Explicit 2x2 multiplication)
    double f00 = F.get(0, 0); double f01 = F.get(0, 1);
    double f10 = F.get(1, 0); double f11 = F.get(1, 1);
    
    double x0_val = x.get(0, 0);
    double x1_val = x.get(1, 0);

    double x0_new = f00 * x0_val + f01 * x1_val;
    double x1_new = f10 * x0_val + f11 * x1_val;

    x.set(0, 0, x0_new);
    x.set(1, 0, x1_new);

    // 2. Covariance prediction
    // temp = P * F^T
    double p00 = P.get(0, 0); double p01 = P.get(0, 1);
    double p10 = P.get(1, 0); double p11 = P.get(1, 1);

    double t00 = p00 * f00 + p01 * f01;
    double t01 = p00 * f10 + p01 * f11;
    double t10 = p10 * f00 + p11 * f01;
    double t11 = p10 * f10 + p11 * f11;

    // P_new = F * temp + Q
    double pn00 = (f00 * t00 + f01 * t10) + Q.get(0, 0);
    double pn01 = (f00 * t01 + f01 * t11) + Q.get(0, 1);
    double pn10 = (f10 * t00 + f11 * t10) + Q.get(1, 0);
    double pn11 = (f10 * t01 + f11 * t11) + Q.get(1, 1);

    P.set(0, 0, pn00); P.set(0, 1, pn01);
    P.set(1, 0, pn10); P.set(1, 1, pn11);
}

void KalmanFilter::update(double z_meas) {
    // Update Step (Correction):
    // K = P * H^T * (H * P * H^T + R)^-1

    // 1. Innovation: y = z - Hx
    // Since H = [0, 1], Hx is simply x[1] (velocity)
    double x_vel = x.get(1, 0);
    double y = z_meas - x_vel;

    // 2. Innovation Covariance: S = H * P * H^T + R
    // For H = [0, 1], H*P*H^T corresponds to P(1,1)
    double p11 = P.get(1, 1);
    double r_val = R.get(0, 0);
    double s = p11 + r_val;

    // 3. Optimal Kalman Gain: K = P * H^T * S^-1
    double k0 = 0.0;
    double k1 = 0.0;
    
    // Numerical stability check
    if (std::abs(s) > 1e-9) {
        k0 = P.get(0, 1) / s;
        k1 = P.get(1, 1) / s;
    }

    // 4. Update State: x = x + K * y
    double x0_new = x.get(0, 0) + k0 * y;
    double x1_new = x.get(1, 0) + k1 * y;
    x.set(0, 0, x0_new);
    x.set(1, 0, x1_new);

    // 5. Update Covariance: P = (I - K * H) * P
    // We compute this explicitly to ensure symmetry and performance
    // (I - KH) matrix elements:
    double m00 = 1.0;    double m01 = -k0;
    double m10 = 0.0;    double m11 = 1.0 - k1;

    double p00 = P.get(0, 0); double p01 = P.get(0, 1);
    double p10 = P.get(1, 0); double p11_old = P.get(1, 1);

    double new_p00 = m00 * p00 + m01 * p10;
    double new_p01 = m00 * p01 + m01 * p11_old;
    double new_p10 = m10 * p00 + m11 * p10;
    double new_p11 = m10 * p01 + m11 * p11_old;

    P.set(0, 0, new_p00); P.set(0, 1, new_p01);
    P.set(1, 0, new_p10); P.set(1, 1, new_p11);
}

double KalmanFilter::getAngle() const { return x.get(0, 0); }
double KalmanFilter::getVelocity() const { return x.get(1, 0); }
Matrix KalmanFilter::getState() const { return x; }