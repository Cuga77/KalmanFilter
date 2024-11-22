#pragma once

#include "matrix.h"
#include "vector.h"

class KalmanFilter {
    Matrix F; // матрица перехода состояния
    Matrix H; // матрица измерений
    Matrix Q; // ковариационная матрица шума системы
    Matrix R; // ковариационная матрица шума измерений
    Matrix P; // ковариационная матрица ошибки оценки
    Matrix K; // усиление Калмана
    my_Vector x; // оценка состояния

public:
    KalmanFilter(double T, double g, double l, double q, double r);

    void initialize(const my_Vector &x0, const Matrix &P0);

    void predict();

    void update(const my_Vector &z);

    my_Vector getState() const;
};
