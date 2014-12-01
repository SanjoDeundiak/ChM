#pragma once

#include <vector>

class TMA1
{
    std::vector<double> alfa, beta;
    std::vector<double> res;
    size_t i; // number of iterations
    size_t n; // size of system
public:
    TMA1(size_t n) : i(0), alfa(n + 1), beta(n + 1), res(n)
    {
        this->n = n;
    }

    std::vector<double> Result();
    void Iterate(double A, double C, double B, double F);
};