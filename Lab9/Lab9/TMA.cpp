#include "TMA.h"

void TMA::Iterate(double A, double B, double C, double D)
{
    if (i == n)
        throw std::exception("Too much iterations");

    if (i == 0)
    {
        c[i] = C / B;
        d[i] = D / B;
    }
    else
    {
        c[i] = C / (B - A*c[i - 1]);
        d[i] = (D - A*d[i - 1]) / (B - A*c[i - 1]);
    }
    i++;
}

std::vector<double> TMA::Result()
{
    if (i != n)
        throw std::exception("Too few iterations");

    r[n - 1] = d[n - 1];;
    for (size_t i = n - 1; i > 0; i--)
        r[i - 1] = d[i - 1] - c[i - 1] * r[i];

    return r;
};