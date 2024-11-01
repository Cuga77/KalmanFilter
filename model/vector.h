// vector.h
#pragma once

#include <vector>
#include <iostream>
#include <cmath>

struct my_Vector {
    std::vector<double> vec;

    my_Vector();
    ~my_Vector();
    my_Vector(int size);
    my_Vector(const my_Vector &copyable);
    my_Vector(my_Vector &&copyable);
    my_Vector &operator =(const my_Vector &copyable);
    my_Vector &operator =(my_Vector &&copyable) noexcept;
    my_Vector reverse();

    friend my_Vector operator +(const my_Vector &v1, const my_Vector &v2);
    friend my_Vector operator -(const my_Vector &v1, const my_Vector &v2);
    friend my_Vector operator *(const my_Vector &v1, double k);
    friend my_Vector operator *(double k, const my_Vector &v1);

    void fill(double a, int count);
    double lorentz_product(my_Vector v2) const;
    double len() const;
    my_Vector cut(int count);
    std::vector<double> get_vec();
    void set(int i, double val);
    void add_coord(int i, double val);
    void print() const;
    void print(std::ostream &os) const;
    double get(int i) const;
};
