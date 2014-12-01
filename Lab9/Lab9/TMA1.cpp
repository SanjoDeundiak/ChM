#include "TMA1.h"

void TMA1::Iterate(double A, double C, double B, double F)
{
    if (i == n)
        throw std::exception("Too much iterations");

    if (i == n - 1 && B != 0)
        throw std::exception("AAA2");

    if (i == 0)
    {
        if (A != 0)
            throw std::exception("AAA3");
        alfa[1] = -B / C;
        beta[1] = F / C;
    }
    else
    {
        alfa[i + 1] = -B / (A*alfa[i] + C);
        beta[i + 1] = (F - A*beta[i]) / (A*alfa[i] + C);
    }
    i++;
    if (i == n)
        res[n - 1] = (F - A*beta[n-1]) / (C + A*alfa[n-1]);
}

std::vector<double> TMA1::Result()
{
    if (i != n)
        throw std::exception("Too few iterations");

    for (size_t i = n - 1; i > 0; i--)
        res[i - 1] = alfa[i] * res[i] + beta[i];

    return res;
};