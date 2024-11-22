#include "vector.h"
#include "matrix.h"

my_Vector::my_Vector() = default;

my_Vector::~my_Vector() {
    vec.clear();
}

my_Vector::my_Vector(int size) {
    vec = std::vector<double>(size);
}

my_Vector::my_Vector(const my_Vector &copyable) {
    vec = copyable.vec;
}

my_Vector::my_Vector(my_Vector &&copyable) {
    vec = std::move(copyable.vec);
}

my_Vector &my_Vector::operator =(const my_Vector &copyable) {
    if (this != &copyable) {
        vec = copyable.vec;
    }
    return *this;
}

my_Vector &my_Vector::operator =(my_Vector &&copyable) noexcept {
    if (this != &copyable) {
        vec = std::move(copyable.vec);
    }
    return *this;
}

my_Vector my_Vector::reverse() {
    my_Vector new_vec(*this);
    std::reverse(new_vec.vec.begin(), new_vec.vec.end());
    return new_vec;
}

my_Vector operator +(const my_Vector &v1, const my_Vector &v2) {
    my_Vector add(v1.vec.size());
    for (size_t i = 0; i < v1.vec.size(); i++) {
        add.vec[i] = v1.vec[i] + v2.vec[i];
    }
    return add;
}

my_Vector operator -(const my_Vector &v1, const my_Vector &v2) {
    my_Vector sub(v1.vec.size());
    for (size_t i = 0; i < v1.vec.size(); i++) {
        sub.vec[i] = v1.vec[i] - v2.vec[i];
    }
    return sub;
}

my_Vector operator *(const my_Vector &v1, double k) {
    my_Vector mul(v1.vec.size());
    for (size_t i = 0; i < v1.vec.size(); i++) {
        mul.vec[i] = v1.vec[i] * k;
    }
    return mul;
}

my_Vector operator *(double k, const my_Vector &v1) {
    return v1 * k;
}

void my_Vector::fill(double a, int count) {
    vec.clear();
    vec.resize(count, a);
}

double my_Vector::lorentz_product(my_Vector v2) const {
    if (vec.size() != v2.vec.size()) return 0;
    double result = 0;
    for (size_t i = 0; i < vec.size() - 1; i++) {
        result += vec[i] * v2.vec[i];
    }
    result -= vec[vec.size() - 1] * v2.vec[v2.vec.size() - 1];
    return result;
}

double my_Vector::len() const {
    double result = 0;
    for (auto i: vec) {
        result += i * i;
    }
    return sqrt(result);
}

my_Vector my_Vector::cut(int count) {
    vec.resize(count);
    return *this;
}

std::vector<double> my_Vector::get_vec() {
    return vec;
}

void my_Vector::set(int i, double val) {
    vec[i] = val;
}

double my_Vector::get(int i) const {
    return vec[i];
}

void my_Vector::add_coord(int i, double val) {
    vec[i] += val;
}

void my_Vector::print() const {
    for (auto i: vec) {
        std::cout << " " << i;
    }
    std::cout << std::endl;
}

void my_Vector::print(std::ostream &os) const {
    for (auto i: vec) {
        os << " " << i;
    }
    os << std::endl;
}
