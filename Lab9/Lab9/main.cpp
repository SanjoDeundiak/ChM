// utt=uxx+F(t,x) 0<x<L=1
// u(t,x) u(0,x)=u0(x), ut(0,x)=0
// u(t,0)=u1(t), u(t,L)=u2(t)
// u0=u(0,x)=(x+0.2)sin(M_PIx/2)
// u1=u(t,0)=0
// u2=u(t,L)=u(t,1)=1.2cos(M_PI*t)
// F(t,x)=-M_PI/2cos(M_PI*t)*(2cos(M_PI*x/2)+M_PI*(x+0.2)sin(M_PI*x/2)
// v0=du/dt(x,0)=0
// ANSWER = (x+0.2)*sin(M_PI*x/2)*cos(M_PI*t)

#include "Hyper.h"
#include "TMA1.h"
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>

func2 answ = [](double t, double x)->double{ return (x + 0.2)*sin(M_PI*x / 2)*cos(M_PI*t); };

int Test()
{
    size_t const Nt = 10;
    size_t const Nl = 5;

    func2 F = [](double t, double x)->double{ return -M_PI*cos(M_PI*t)*(cos(M_PI*x / 2) + 3 / 4 * M_PI*(x + 0.2)*sin(M_PI*x / 2)); };
    func u0 = [](double x)->double{ return (x + 0.2)*sin(M_PI*x / 2); };
    func u1 = [](double t)->double{ return 0.; };
    func u2 = [](double t)->double{ return 1.2*cos(M_PI*t); };
    func v0 = [](double x)->double{ return 0.; };

    double period = 2, length = 1;
    double h = length / Nl;
    double dt = period / Nt;

    std::vector<std::vector<double>> res;
    try
    {
        res = Hyper::Process(0.5, F, u0, u1, u2, v0, Nt, Nl, period, length);
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        std::cin.get();
        return -1;
    }

    std::ofstream out("output.txt");

    TMA1 tma(3);
    tma.Iterate(0, 1, 0, 1);
    tma.Iterate(1, 1, 1, 6);
    tma.Iterate(2, 1, 0, 7);
    std::vector<double> r = tma.Result();
    for (double a : r)
        out << a << ' ';
    out << std::endl;

    // TestMax
    out << std::left;
    {
        for (size_t k = 7; k == 7; k++)
        {
            for (size_t n = 0; n < Nl + 1; n++)
            {
                out << std::setw(4) << n*h << '\t';
                out << std::setw(12) << answ(k*dt, n*h) << '\t' << std::setw(12) << res[k][n] << std::endl;//res[k][n] << '\t';
            }

            out << std::endl;
        }
    }

    out << std::endl;

    std::cin.get();

    return 0;
}

int main()
{
    return Test();
}