#include "RKM.h"

std::function<double(double)> RKM::Solution(const std::function<double(double, double)>& func, double x0, double y0, int n /*=10*/)
{
    RKMResult result;
    result.h = 1./abs(n);
    result.x0 = x0;
    result.y0 = y0;
    result.deriv = func;

    return result;
}

double RKM::RKMResult::operator()(double xr)
{
    if (xr > x0)
        h = abs(h);
    else
        h = -abs(h);

    double k1, k2, k3, k4;
    double x = x0, y = y0;

    while (abs(x - xr) > abs(h / 2))
    {
        k1 = deriv(x, y);
        k2 = deriv(x + h / 2, y + k1 * h / 2);
        k3 = deriv(x + h / 2, y + k2 * h / 2);
        k4 = deriv(x + h, y + k3 * h);

        y += (k1 + 2*k2 + 2*k3 + k4) * h/6;
        x += h;
    }

    return y;
}

std::vector<double> RKM::Result(const std::function<double(double, double)>& deriv, double x0, double y0, int n /*= 10*/)
{
    std::vector<double> result{ y0 };
    double h = 1./abs(n);

    double k1, k2, k3, k4;
    double x = x0, y = y0;

    while (abs(x - x0) < 1)
    {
        k1 = deriv(x, y);
        k2 = deriv(x + h / 2, y + k1 * h / 2);
        k3 = deriv(x + h / 2, y + k2 * h / 2);
        k4 = deriv(x + h, y + k3 * h);

        y += (k1 + 2 * k2 + 2 * k3 + k4) * h / 6;
        x += h;
        result.push_back(y);
    }

    return result;
}
