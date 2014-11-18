#pragma once

#include <vector>

class TMA
{
    std::vector<double> c, d, r;
    size_t i; // number of iterations
    size_t n; // size of system
public:
    TMA(size_t n) : i(0), c(n), d(n), r(n)
    {
        this->n = n;
    }

    std::vector<double> Result();
    void Iterate(double a, double b, double c, double d);
};