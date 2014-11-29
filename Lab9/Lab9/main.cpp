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
#include <iostream>
#include <algorithm>
#include <fstream>

size_t const Nt = 10;
size_t const Nl = 20;
double const c = 1.;

int main()
{
    func2 F = [](double t, double x)->double{ return -M_PI*cos(M_PI*t)*(cos(M_PI*x/2) + 3/4*M_PI*(x + 0.2)*sin(M_PI*x / 2)); };
    func u0 = [](double x)->double{ return (x + 0.2)*sin(M_PI*x/2); };
    func u1 = [](double t)->double{ return 0.; };
    func u2 = [](double t)->double{ return 1.2*cos(M_PI*t); };
    func v0 = [](double x)->double{ return 0.; };
    func2 answ = [](double t, double x)->double{ return (x + 0.2)*sin(M_PI*x / 2)*cos(M_PI*t); };

    std::vector<std::vector<double>> res;
    try
    {
        res = Hyper::Process(c, 0.5, F, u0, u1, u2, v0, Nt, Nl);
    }
    catch (std::exception& e)
    {
        std::cout << e.what() << std::endl;
        std::cin.get();
        return -1;
    }

    std::ofstream out("output.txt");

    for (size_t k = 0; k < Nt + 1; k++)
    {
        size_t i;
        double m = 0;
        for (size_t n = 0; n < Nl + 1; n++)
        {
            double z = abs(res[k][n] - answ(2.*k / Nt, 1.*n / Nl));
            if (z > m)
            {
                m = z;
                i = n;
            }
        }
        out << k*2. / Nt << '\t' << m << '\t' << i << std::endl;
        //out << answ(2.*k/Nt, 1) << '\t' << res[k][Nl] << '\t' << abs(answ(2.*k / Nt, 1) - res[k][Nl]) << std::endl;
    }
    
    /*for (size_t n = 0; n < Nl; n++)
    {
        out << answ(2. / Nt, 1.*n/Nl) << '\t' << res[1][n] << '\t' << abs(answ(2. / Nt, 1.*n/Nl) - res[1][n]) << std::endl;
    }*/

    out << std::endl;

    std::cin.get();
}