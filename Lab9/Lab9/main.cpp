#include "Hyper.h"
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>

int Test()
{
    size_t const Nt = 1000;
    size_t const Nl = 1000;

    func u0 = [](double x)->double{ return (x + 0.2)*sin(M_PI*x / 2); };
    func u1 = [](double t)->double{ return 0.; };
    func u2 = [](double t)->double{ return 1.2*cos(M_PI*t); };
    func v0 = [](double x)->double{ return 0.; };

    func2 utt = [](double t, double x)->double{ return -M_PI*M_PI*(x + 0.2)*sin(M_PI*x / 2)*cos(M_PI*t); };
    func2 uxx = [](double t, double x)->double{ return (-M_PI*M_PI/4*(x+0.2)*sin(M_PI*x/2)+M_PI*cos(M_PI*x/2))*cos(M_PI*t); };

    func2 F = [utt, uxx](double t, double x)->double{ return utt(t,x)-uxx(t,x); };

    func2 answ = [](double t, double x)->double{ return (x + 0.2)*sin(M_PI*x / 2)*cos(M_PI*t); };

    double period = 2, length = 1;
    double h = length / Nl;
    double dt = period / Nt;

    std::vector<std::vector<double>> res = Hyper::Process(0.75, F, u0, u1, u2, v0, Nt, Nl, period, length);

    std::ofstream out("output.txt");

    out << std::left << std::fixed;

    {
        for (size_t k = 0; k < Nt + 1; k++)
        {
            double m = 0;
            for (size_t n = 0; n < Nl + 1; n++)
            {
                double a1 = answ(k*dt, n*h);
                double a2 = res[k][n];
                double err = abs(a1-a2);
                m = std::max(m,err);
            }
            out << std::setw(4) << k*dt << '\t'
                << std::setw(12) << m << std::endl;
        }
    }

    /*{
        for (size_t k = 0; k < Nt + 1; k++)
        {
            double m = 0;
            for (size_t n = 0; n < Nl + 1; n++)
            if (abs(res[k][n]) > m)
                m = abs(res[k][n]);

            out << std::setw(4) << k*dt << '\t'
                << std::setw(12) << m << std::endl;
        }
    }*/

    /*{
        for (size_t n = 0; n < Nl + 1; n++)
        out << std::setw(4) << n*h << '\t'
            << std::setw(12) << res[Nt][n] << std::endl;
    }*/

    /*{
        for (size_t n = 0; n < Nl + 1; n++)
            out << std::setw(4) << n*h << '\t'
            << std::setw(12) << answ(Nt*dt, n*h) - res[Nt][n] << std::endl;
    }*/

    /*for (size_t k = 0; k <= Nt; k++)
    {
        out << "t=" << k*dt << std::endl;
        for (size_t n = 0; n <= Nl; n++)
        {
            out << n*h << '\t' << answ(k*dt, n*h) << '\t' << res[k][n] << '\t' << answ(k*dt, n*h) - res[k][n] << std::endl;
        }
        out << std::endl;
    }*/

    out << std::endl;

    return 0;
}

int main()
{
    return Test();
}