#include "Hyper.h"
#include "TMA.h"

extern func2 answ;

#include <fstream>
#include <algorithm>
std::ofstream sf("s.txt");

std::vector<std::vector<double>> Hyper::Process(double s, func2 F,
    func u0, func u1, func u2, func v0,
    size_t Nt, size_t Nl, double period, double length)
{
    std::vector<std::vector<double>> u(Nt+1, std::vector<double>(Nl+1)); // first time, next x

    double h = length / Nl;
    double dt = period / Nt;

    double q = h / dt*h / dt;

    for (size_t n = 0; n <= Nl; n++)
        u[0][n] = u0(n*h);

    u[1][0] = u1(dt);
    u[1][Nl] = u2(dt);
    for (size_t n = 1; n < Nl; n++)
        u[1][n] = u[0][n] + dt*v0(h*n) + dt*dt / 2 * ((u[0][n + 1] - 2 * u[0][n] + u[0][n - 1]) / h / h +F(0, h*n));
        //u[1][n] = answ(dt, n*h);

    double err = 0;
    for (size_t k = 1; k < Nt; k++) // time loop
    {
        TMA tma(Nl + 1);
        tma.Iterate(0, 1, 0, u1((k + 1)*dt));

        for (size_t n = 1; n < Nl; n++)
        {
            tma.Iterate(-s, q + 2 * s, -s,
                2 * u[k][n]*q - u[k - 1][n]*q + (1 - 2 * s)*(u[k][n + 1] - 2 * u[k][n] + u[k][n - 1]) + s*(u[k - 1][n + 1] - 2 * u[k - 1][n] + u[k - 1][n - 1])
                + (s*F((k + 1)*dt, n*h) + (1 - 2 * s)*F(k*dt, n*h) + s*F((k - 1)*dt, n*h))*h*h);
                //2 * answ(k*dt, n*h) * q - answ((k - 1)*dt, n*h) * q + (1 - 2 * s)*(answ(k*dt, (n + 1)*h) - 2 * answ(k*dt, n*h) + answ(k*dt, (n - 1)*h))
                //+ s*(answ((k - 1)*dt, (n + 1)*h) - 2 * answ((k - 1)*dt, n*h) + answ((k - 1)*dt, (n - 1)*h))
                //+ F(k*dt, n*h) *h*h);
        }
        tma.Iterate(0, 1, 0, u2((k + 1)*dt));

        u[k+1] = tma.Result();
    }

    return u;
}