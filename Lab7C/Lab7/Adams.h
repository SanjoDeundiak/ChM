#pragma once

#include <functional>
#include <vector>

class Adams
{
public:
    static std::vector<double> Result(const std::function<double(double, double)>& deriv, double x0, double y0, int n = 10);
};