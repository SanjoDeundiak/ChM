#include "Adams.h"
#include "RKM.h"

std::vector<double> Adams::Result(const std::function<double(double, double)>& deriv, double x0, double y0, int n /*= 10*/)
{
    double h = 1./abs(n);
    double f[4];
    std::vector<double> y;

    // Getting initial values using Runge-Kutta method
    double x = x0;
    std::function<double(double)> rk = RKM::Solution(deriv, x0, y0, n);
    for (size_t i = 0; i < 4; i++)
    {
        y.push_back(rk(x));
        x += h;
    }

    for (int i = 4; i <= n; i++)
    {
        y.push_back(y[i-1] + h/24*(55*deriv(x - h, y[i-1]) - 59*deriv(x - 2*h, y[i-2]) + 37*deriv(x - 3*h, y[i-3]) - 9*deriv(x - 4*h, y[i-4])));
        x += h;
    }

    return y;
}
