#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include "matrix.h"

Matrix::Matrix(int n_new, int m_new) {
    n = n_new;
    m = m_new;
    for (int i = 0; i < n; i++) {
        std::vector<double> vec(m);
        matrix.push_back(vec);
    }
}

Matrix::Matrix(int n_new, int m_new, double val) {
    n = n_new;
    m = m_new;
    for (int i = 0; i < n; i++) {
        std::vector<double> vec(m);
        vec[i] = val;
        matrix.push_back(vec);
    }
}

Matrix::Matrix(const std::string &file) {
    std::ifstream ifs;
    ifs.open(file);
    ifs >> m >> n;
    for (int i = 0; i < n; i++) {
        std::vector<double> row;
        for (int j = 0; j < m; j++) {
            double elem;
            ifs >> elem;
            row.push_back(elem);
        }
        matrix.push_back(row);
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < n; i++) {
        matrix[i].clear();
    }
    matrix.clear();
    n = 0;
    m = 0;
}

Matrix::Matrix(const Matrix &copy) {
    n = copy.n;
    m = copy.m;
    for (int i = 0; i < n; i++) {
        std::vector<double> v;
        for (int j = 0; j < m; j++) {
            v.push_back(copy.matrix[i][j]);
        }
        matrix.push_back(v);
    }
}

Matrix::Matrix(Matrix &&copy) {
    n = copy.n;
    m = copy.m;
    for (int i = 0; i < n; i++) {
        std::vector<double> v;
        for (int j = 0; j < m; j++) {
            v.push_back(copy.matrix[i][j]);
        }
        matrix.push_back(v);
    }
}

Matrix &Matrix::operator=(const Matrix &copy) {
    for (int i = 0; i < n; i++) {
        matrix[i].clear();
    }
    matrix.clear();
    m = copy.m;
    n = copy.n;
    for (int i = 0; i < n; i++) {
        std::vector<double> v;
        for (int j = 0; j < m; j++) {
            v.push_back(copy.matrix[i][j]);
        }
        matrix.push_back(v);
    }
    return *this;
}

Matrix &Matrix::operator=(Matrix &&copy) {
    for (int i = 0; i < n; i++) {
        matrix[i].clear();
    }
    matrix.clear();
    m = copy.m;
    n = copy.n;
    for (int i = 0; i < n; i++) {
        std::vector<double> v;
        for (int j = 0; j < m; j++) {
            v.push_back(copy.matrix[i][j]);
        }
        matrix.push_back(v);
    }
    return *this;
}

Matrix operator*(const Matrix &A, double k) {
    Matrix m(A);
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            m.matrix[i][j] *= k;
        }
    }
    return m;
}

Matrix operator*(double k, const Matrix &A) {
    Matrix m(A);
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            m.matrix[i][j] *= k;
        }
    }
    return m;
}

Matrix operator+(const Matrix &B, const Matrix &A) {
    Matrix m(A);
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            m.matrix[i][j] += B.matrix[i][j];
        }
    }
    return m;
}

Matrix operator-(const Matrix &B, const Matrix &A) {
    Matrix m(A);
    for (int i = 0; i < m.n; i++) {
        for (int j = 0; j < m.m; j++) {
            m.matrix[i][j] -= B.matrix[i][j];
        }
    }
    return m;
}

my_Vector operator*(const Matrix &A, const my_Vector &v) {
    my_Vector res(A.n);
    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < A.m; j++) {
            res.vec[i] += A.matrix[i][j] * v.vec[j];
        }
    }
    return res;
}

my_Vector operator*(const my_Vector &v, const Matrix &A) {
    my_Vector res(A.n);
    res.fill(0, v.vec.size());
    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < A.m; j++) {
            res.vec[i] += A.matrix[j][i] * v.vec[j];
        }
    }
    return res;
}

Matrix operator*(const Matrix &A, const Matrix &B) {
    Matrix res(A.n, B.m);
    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < B.m; j++) {
            double sum = 0;
            for (int k = 0; k < A.m; k++) {
                sum += A.matrix[i][k] * B.matrix[k][j];
            }
            res.matrix[i][j] = sum;
        }
    }
    return res;
}

Matrix Matrix::cut(int row, int column) const {
    Matrix new_matrix(n - 1, m - 1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (i != row && j != column) {
                int r = i > row ? i - 1 : i;
                int c = j > column ? j - 1 : j;
                new_matrix.matrix[r][c] = matrix[i][j];
            }
        }
    }
    return new_matrix;
}

Matrix Matrix::transpose() const {
    Matrix new_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            new_matrix.matrix[i][j] = matrix[j][i];
        }
    }
    return new_matrix;
}

double Matrix::determinant() const {
    if (n == m && n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double result = 0;
    for (int i = 0; i < m; i++) {
        result += matrix[0][i] * (i % 2 == 0 ? 1 : -1) * (this->cut(0, i)).determinant();
    }
    return result;
}

Matrix Matrix::choleskyDecomposition() const {
    Matrix L(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            if (j == i) {
                for (int k = 0; k < j; k++) {
                    sum += L.matrix[j][k] * L.matrix[j][k];
                }
                L.matrix[j][j] = sqrt(matrix[j][j] - sum);
            } else {
                for (int k = 0; k < j; k++) {
                    sum += L.matrix[i][k] * L.matrix[j][k];
                }
                L.matrix[i][j] = (matrix[i][j] - sum) / L.matrix[j][j];
            }
        }
    }
    return L;
}

my_Vector Matrix::solve_system(const my_Vector &b) const {
    my_Vector x(n);
    Matrix L = this->choleskyDecomposition();

    // L*y = b
    my_Vector y(n);
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += y.vec[j] * L.matrix[i][j];
        }
        y.vec[i] = (b.vec[i] - sum) / L.matrix[i][i];
    }

    // L^T*x = y
    Matrix Lt = L.transpose();
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += x.vec[j] * Lt.matrix[i][j];
        }
        x.vec[i] = (y.vec[i] - sum) / Lt.matrix[i][i];
    }

    return x;
}

void Matrix::set(int i, int j, double val) {
    matrix[i][j] = val;
}

double Matrix::get(int i, int j) const {
    return matrix[i][j];
}

Matrix Matrix::inverseDiagonal() const {
    Matrix new_matrix(*this);
    for (int i = 0; i < this->matrix.size(); i++) {
        if (this->get(i, i) != 0)
            new_matrix.set(i, i, 1 / this->get(i, i));
    }
    return new_matrix;
}

double Matrix::get_minor(int i, int j) const {
    return cut(i, j).determinant();
}

double Matrix::get_alg(int i, int j) const {
    return ((i + j) % 2 == 0 ? 1 : -1) * get_minor(i, j);
}

Matrix Matrix::invert() const {
    if (determinant() == 0) std::cout << "det = 0" << std::endl;
    Matrix alg(n, m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            alg.matrix[i][j] = get_alg(i, j);
        }
    }
    return 1.0 / determinant() * (alg.transpose());
}

void Matrix::swap(int s1, int s2) {
    std::vector tmp(matrix[s1]);
    matrix[s1] = matrix[s2];
    matrix[s2] = tmp;
}

int Matrix::get_max_index_in_column(int column) const {
    int res = -1;
    double max = 0;
    for (int i = column; i < n; i++) {
        if (fabs(matrix[i][column]) > max) {
            max = fabs(matrix[i][column]);
            res = i;
        }
    }
    return res;
}

void Matrix::mod_str(double k, int from, int to) {
    for (int i = 0; i < m; i++) {
        matrix[to][i] += matrix[from][i] * k;
    }
}

void Matrix::get_upper_triangle() {
    for (int col = 0; col < m; col++) {
        int mx = get_max_index_in_column(col);
        if (mx != -1) {
            swap(col, mx);
            for (int row = col + 1; row < n; row++) {
                mod_str(-matrix[row][col] / matrix[col][col], col, row);
            }
        }
    }
}

void Matrix::print() const {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

Matrix GetIdentity(int n, int m) {
    Matrix matrix(n, m, 0);
    for (int i = 0; i < n; i++)
        matrix.set(i, i, 1);
    return matrix;
}

Matrix getQMatrix(double dt, double sigma, Matrix lastMtx) {
    double p = dt * dt * dt * dt / 4;
    double f = dt * dt * dt / 2;
    double a = dt * dt / 2;
    double b = a * 2;

    double mtxPattern[9][9] = {
        {p, 0, 0, f, 0, 0, a, 0, 0},
        {0, p, 0, 0, f, 0, 0, a, 0},
        {0, 0, p, 0, 0, f, 0, 0, a},
        {f, 0, 0, b, 0, 0, dt, 0, 0},
        {0, f, 0, 0, b, 0, 0, dt, 0},
        {0, 0, f, 0, 0, b, 0, 0, dt},
        {a, 0, 0, dt, 0, 0, 1, 0, 0},
        {0, a, 0, 0, dt, 0, 0, 1, 0},
        {0, 0, a, 0, 0, dt, 0, 0, 1},
    };

    for (int i = 0; i < lastMtx.matrix.size(); i++) {
        for (int j = 0; j < lastMtx.matrix[i].size(); j++) {
            lastMtx.set(i, j, mtxPattern[i][j]);
        }
    }
    return lastMtx * (sigma * sigma);
}

Matrix getRMatrix(double sigma) {
    return GetIdentity(3, 3) * (sigma * sigma);
}

Matrix getPMatrix(double sigma) {
    return GetIdentity(9, 9) * (sigma * sigma);
}

Matrix getFMatrix(double t, Matrix lastMtx) {
    double a = t * t / 2;
    double mtxPettern[9][9] = {
        {1, 0, 0, t, 0, 0, a, 0, 0},
        {0, 1, 0, 0, t, 0, 0, a, 0},
        {0, 0, 1, 0, 0, t, 0, 0, a},
        {0, 0, 0, 1, 0, 0, t, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, t, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, t},
        {0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1},
    };

    for (int i = 0; i < lastMtx.matrix.size(); i++) {
        for (int j = 0; j < lastMtx.matrix[i].size(); j++) {
            lastMtx.set(i, j, mtxPettern[i][j]);
        }
    }
    return lastMtx;
}

Matrix getHMatrix() {
    Matrix m(3, 9, 0);
    for (int i = 0; i < 3; i++)
        m.set(i, i, 1);
    return m;
}
