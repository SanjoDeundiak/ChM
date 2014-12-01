#include "Hyper.h"
#include "TMA1.h"

extern func2 answ;

#include <fstream>
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
        //u[1][n] = u[0][n] + dt*v0(h*n) + dt*dt / 2 * ((u[0][n + 1] - 2 * u[0][n] + u[0][n - 1]) / h / h +F(0, h*n));
        u[1][n] = answ(dt, n*h);

    for (size_t k = 1; k < Nt; k++) // time loop
    {
        TMA1 tma(Nl + 1);
        std::vector<std::vector<double>> A(Nl+1, std::vector<double>(Nl+1));
        std::vector<double> B(Nl+1);
        tma.Iterate(0, 1, 0, u1((k + 1)*dt));
        A[0][0] = 1;
        B[0] = u1((k + 1)*dt);
        for (size_t i = 1; i <= Nl; i++)
            A[0][i] = 0;

        for (size_t n = 1; n < Nl; n++)
        {
            for (size_t i = 0; i < n - 1; i++)
                A[n][i] = 0;
            for (size_t i = n + 2; i < Nl; i++)
                A[n][i] = 0;

            A[n][n - 1] = -s;
            A[n][n] = q + 2*s;
            A[n][n + 1] = -s;

            tma.Iterate(-s, q + 2 * s, -s,
                2 * u[k][n]*q - u[k - 1][n]*q + (1 - 2 * s)*(u[k][n + 1] - 2 * u[k][n] + u[k][n - 1]) + s*(u[k - 1][n + 1] - 2 * u[k - 1][n] + u[k - 1][n - 1])
                + (s*F((k + 1)*dt, n*h) + (1 - 2 * s)*F(k*dt, n*h) + s*F((k - 1)*dt, n*h))*h*h);
            
            B[n] = 2 * u[k][n] * q - u[k - 1][n] * q + (1 - 2 * s)*(u[k][n + 1] - 2 * u[k][n] + u[k][n - 1]) + s*(u[k - 1][n + 1] - 2 * u[k - 1][n] + u[k - 1][n - 1])
                + (s*F((k + 1)*dt, n*h) + (1 - 2 * s)*F(k*dt, n*h) + s*F((k - 1)*dt, n*h))*h*h;
        }
        tma.Iterate(0, 1, 0, u2((k + 1)*dt));
        for (size_t i = 0; i < Nl; i++)
            A[Nl][i] = 0;
        A[Nl][Nl] = 1;
        B[Nl] = u2((k + 1)*dt);

        if (k == 7)
        {
            for (size_t i = 0; i < Nl+1; i++)
            {
                sf << '{';
                for (size_t j = 0; j < Nl+1; j++)
                {
                    sf << A[i][j] << ", ";
                }
                sf << '}';
            }
            sf << std::endl;
            sf << '{';
            for (size_t i = 0; i < Nl + 1; i++)
            {
                sf << B[i] << ", ";
            }
            sf << '}';
        }

        u[k+1] = tma.Result();
        sf << std::endl;
        for (size_t i = 0; k == 7 && i <= Nl; i++)
            sf << u[k + 1][i] << '\t';
    }

    return u;
}