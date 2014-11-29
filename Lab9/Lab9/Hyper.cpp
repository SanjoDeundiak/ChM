#include "Hyper.h"
#include "TMA.h"

std::vector<std::vector<double>> Hyper::Process(double c, double s, func2 F,
    func u0, func u1, func u2, func v0,
    size_t Nt, size_t Nl)
{
    std::vector<std::vector<double>> u(Nt+1, std::vector<double>(Nl+1)); // first time, next x

    double h = 1 / Nl;
    double dt = 2 / Nt;

    double q = dt*dt / h / h;

    for (size_t n = 0; n <= Nl; n++)
        u[0][n] = u0(n*h);

    u[1][0] = u1(dt);
    u[1][Nl] = u2(dt);
    for (size_t n = 1; n < Nl; n++)
        u[1][n] = u[0][n] + dt*v0(h*n) + dt*dt / 2 * (c*c*(u[0][n+1]-2*u[0][n]+u[0][n-1]) / h / h + F(0, h*n));

    for (size_t k = 1; k < Nt; k++) // time loop
    {
        TMA tma(Nt - 1);
        u[k][0] = u1(k*dt);
        u[k][Nl] = u2(k*dt);

        for (size_t n = 2; n < Nt; n++)
        {
            tma.Iterate(-c*c*s*q, 1 + 2 * c*c*q, -c*c*s*q,
                2 * u[k][n] - u[k - 1][n] + c*c*(1 - 2 * s)*q*(u[k][n + 1] - 2 * u[k][n] + u[k][n - 1]) + c*c*s*q*(u[k - 1][n + 1] - 2 * u[k - 1][n] + u[k - 1][n - 1])
                + dt*dt*(s*F((k + 1)*dt, (n + 1)*h) + (1 - 2 * s)*F(k*dt, n*h) + s*F((k - 1)*dt, n*h)));
        }
        std::vector<double> v = tma.Result();
        u[k].insert(u[k].begin() + 1, v.begin(), v.end());
    }

    return u;
}