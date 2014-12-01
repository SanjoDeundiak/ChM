#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <functional>
#include <vector>

typedef std::function<double(double)> func;
typedef std::function<double(double, double)> func2;

class Hyper
{
public:
    static std::vector<std::vector<double>> Process(double s, func2 F, func u0, func u1, func u2, func v0, size_t Nt, size_t Nl, double period, double length);
};