// #include "KalmanFilter.h"


//todo: все очент грустно


//
// KalmanFilter::KalmanFilter(double T, double g, double l, double q, double r) {
//     F = Matrix(2, 2);
//     F.set(0, 0, 1);
//     F.set(0, 1, T);
//     F.set(1, 0, - (g / l) * T);
//     F.set(1, 1, 1);
//
//     H = Matrix(1, 2);
//     H.set(0, 0, 1);
//     H.set(0, 1, 0);
//
//     Q = Matrix(2, 2);
//     Q.set(0, 0, q * (T*T*T)/3.0);
//     Q.set(0, 1, q * (T*T)/2.0);
//     Q.set(1, 0, q * (T*T)/2.0);
//     Q.set(1, 1, q * T);
//
//     R = Matrix(1, 1);
//     R.set(0, 0, r);
// }
//
// void KalmanFilter::initialize(const my_Vector& x0, const Matrix& P0) {
//     x = x0;
//     P = P0;
// }
//
// void KalmanFilter::predict() {
//     x = F * x;
//     P = F * P * F.transpose() + Q;
// }
//
// void KalmanFilter::update(const my_Vector& z) {
//     const Matrix S = H * P * H.transpose() + R;
//     K = P * H.transpose() * S.invert();
//
//     my_Vector y = z - (H * x);
//     x = x + (K * y);
//
//     Matrix I = GetIdentity(P.n, P.m);
//     P = (I - (K * H)) * P;
// }
//
// my_Vector KalmanFilter::getState() const {
//     return x;
// }
