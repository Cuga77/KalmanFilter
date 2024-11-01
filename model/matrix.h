// matrix.h
#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include "vector.h"

struct Matrix {
    std::vector<std::vector<double>> matrix;
    int n = 0, m = 0;

    Matrix();
    Matrix(int n_new, int m_new);
    Matrix(int n_new, int m_new, double val);
    Matrix(const std::string &file);
    ~Matrix();
    Matrix(const Matrix &copy);
    Matrix(Matrix &&copy);
    Matrix &operator =(const Matrix &copy);
    Matrix &operator =(Matrix &&copy);

    friend Matrix operator *(const Matrix &A, double k);
    friend Matrix operator *(double k, const Matrix &A);
    friend Matrix operator +(const Matrix &B, const Matrix &A);
    friend Matrix operator -(const Matrix &B, const Matrix &A);
    friend my_Vector operator *(const Matrix &A, const my_Vector &v);
    friend my_Vector operator *(const my_Vector &v, const Matrix &A);
    friend Matrix operator *(const Matrix &A, const Matrix &B);

    Matrix cut(int row, int column) const;
    double determinant() const;
    void set(int i, int j, double val);
    double get(int i, int j) const;
    Matrix inverseDiagonal() const;
    double get_minor(int i, int j) const;
    double get_alg(int i, int j) const;
    Matrix transpose() const;
    Matrix invert() const;
    void swap(int s1, int s2);
    int get_max_index_in_column(int column) const;
    void mod_str(double k, int from, int to);
    void get_upper_triangle();
    Matrix choleskyDecomposition() const;
    my_Vector solve_system(my_Vector b) const;
    void print();
};

Matrix GetIdentity(int n, int m);
Matrix getQMatrix(double dt, double sigma, Matrix lastMtx);
Matrix getRMatrix(double sigma);
Matrix getPMatrix(double sigma);
Matrix getFMatrix(double t, Matrix lastMtx);
Matrix getHMatrix();
