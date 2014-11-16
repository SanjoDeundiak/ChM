#pragma once

#include <vector>

std::vector<double> operator+(const std::vector<double> &v1, const std::vector<double> &v2);
std::vector<double> operator*(const std::vector<double> v, double n);
std::vector<double> operator*(double n, const std::vector<double> v);
std::ostream &operator<<(std::ostream &os, const std::vector<double> &v);

template <typename T>
std::vector<T>& reverse(std::vector<T>& v)
{
    for (size_t i = 0; i < v.size() / 2; i++)
        std::swap(v[i], v[v.size() - i - 1]);

    return v;
}