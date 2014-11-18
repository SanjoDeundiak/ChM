#pragma once

#include <vector>
#include <functional>

typedef std::function<double(double)> func;

class FDM
{
public:
    static std::vector<double> Result(func p, func q, func f,
        double a, double kda, double ka, double f1,
        double b, double kdb, double kb, double f2,
        size_t N);
};