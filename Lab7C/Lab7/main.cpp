// y' = (2x + y - 3) / (x - 1)
// y(2) = 1
// answer: y = 2(x-1)log(x-1) + 1

#include "RKM.h"
#include "Adams.h"
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>

int main()
{
    std::function<double(double)> answer = [](double x)->double{return 2 * (x - 1)*log(x - 1) + 1; };
    std::function<double(double, double)> deriv = [](double x, double y)->double{ return (2 * x + y - 3) / (x - 1); };

    std::ofstream output("output.txt");

    // Runge-Kutta
    {
        int n = 10;
        double x = 2, y = 1;
        std::vector<double> res = RKM::Result(deriv, 2, 1);

        double h = 1. / abs(n);
        output << "RUNGE-KUTTA h = " << h << std::endl << std::endl;
        for (double y : res)
        {
            output << x << '\t' << y << '\t' << answer(x) << '\t' << abs(y - answer(x)) << std::endl;
            x += h;
        }
    }

    output << std::endl;

    // Adams
    {
        int n = 10;
        double x = 2, y = 1;
        std::vector<double> res = Adams::Result(deriv, x, y, n);

        double h = 1./abs(n);
        output << "ADAMS h = " << h << std::endl << std::endl;
        for (double y : res)
        {
            output << x << '\t' << y << '\t' << answer(x) << '\t' << abs(y - answer(x)) << std::endl;
            x += h;
        }
    }
}