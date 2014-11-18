#include "FDM.h"
#include "TMA.h"

#include <iostream>

std::vector<double> FDM::Result(func p, func q, func f,
    double a, double kda, double ka, double f1,
    double b, double kdb, double kb, double f2,
    size_t N)
{
    TMA tma(N + 1);
    double h = (b - a) / N, x = a + h;

    tma.Iterate(0., -kda / h + ka, kda / h, f1);

    for (size_t i = 1; i < N; i++, x += h)
        tma.Iterate(1 / (h*h) - p(x) / (2 * h), -2 / (h*h) + q(x), 1 / (h*h) + p(x) / (2 * h), f(x));

    tma.Iterate(-kdb/h, kdb/h + kb, 0, f2);
    return tma.Result();
}
