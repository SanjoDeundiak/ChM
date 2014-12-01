// Ay = y'' - y' + 2y / x
// a1y = y(1, 1) - 0.5y'(1, 1)
// a2y = y'(1, 4)
// -0.282, 1.981, 2.065, -0.577, -2.029
// y(x) = ax^2 +bx + c + (dx+e)^-1

#include <functional>
#include <fstream>
#include <iostream>

#include "FDM.h"

using namespace std;

const double a = -0.282, b = 1.981, c = 2.065, d = -0.577, e = -2.029;

void test()
{
    func y = [](double x)->double {return a*x*x + b*x + c + 1 / (d*x + e); };
    func yd = [](double x)->double {return 2 * a*x + b - d / pow(d*x + e, 2); };

    func f = [](double x)->double {return 2 * a - b
        + 2 * d*d / pow(d*x + e, 3)
        - 2 * a*x + d / pow(d*x + e, 2)
        + 2 / x*(a*x*x + b*x + c + 1 / (d*x + e)); };
    func p = [](double x)->double {return -1; }, q = [](double x)->double {return 2 / x; };

    double kda = -0.5, ka = 1;
    double kdb = 1, kb = 0;

    double xa = 1.1, f1 = ka*y(xa) + kda*yd(xa);
    double xb = 1.4, f2 = kb*y(xb) + kdb*yd(xb);

    ofstream err("error.txt");
    for (size_t N = 9; N == 9; N *= 2)
    {
        vector<double> res = FDM::Result(p, q, f, xa, kda, ka, f1, xb, kdb, kb, f2, N);
        double x = xa, h = (xb - xa) / N;
        double error = abs(y(x) - res[0]);
        for (double yr : res)
        {
            double t = abs(y(x) - yr);
            if (t > error)
                error = t;
            x += h;
        }
        err << h << '\t' << error << endl;
    }
}

int main()
{
    test();

    cin.get();
}