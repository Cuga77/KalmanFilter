#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include "matrix.h"

class KalmanFilter {
private:
    Matrix x; // State vector (2x1)
    Matrix P; // Covariance matrix (2x2)
    Matrix F; // State Transition matrix (2x2)
    Matrix Q; // Process Noise covariance (2x2)
    Matrix R; // Measurement Noise covariance (1x1)
    Matrix H; // Observation matrix (1x2)

public:
    // Constructor
    KalmanFilter(Matrix x0, Matrix P0, Matrix F, Matrix Q, Matrix R, Matrix H);

    // Setters
    void setF(const Matrix& newF);

    // Prediction Step: x = Fx, P = FPF^T + Q
    void predict();

    // Update Step: x = x + K(z - Hx)
    void update(double z_measurement);

    // Getters
    double getAngle() const;
    double getVelocity() const;
    Matrix getState() const;
};

#endif //KALMANFILTER_H