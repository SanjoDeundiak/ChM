#pragma once

#include <functional>
#include <vector>

class RKM
{
    struct RKMResult
    {
        double h;
        std::function<double(double, double)> deriv;
        double x0, y0;

        double operator()(double x);
    };

public:
    static std::function<double(double)> Solution(const std::function<double(double, double)>& func, double x0, double y0, int n = 10);
    static std::vector<double> Result(const std::function<double(double, double)>& deriv, double x0, double y0, int n = 10);
};