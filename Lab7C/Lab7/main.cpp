// y' = (2x + y - 3) / (x - 1)
// y(2) = 1
// answer: y = 2(x-1)log(x-1) + 1

#include "RKM.h"
#include "Adams.h"
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <iomanip>

int main()
{
    std::function<double(double)> answer = [](double x)->double{return 2 * (x - 1)*log(x - 1) + 1; };
    std::function<double(double, double)> deriv = [](double x, double y)->double{ return (2 * x + y - 3) / (x - 1); };

    std::ofstream output("output.txt");

    for (int n = 10; n <= 100; n += 10)
    {
        double x = 2, y = 1;
        std::vector<double> res1 = RKM::Result(deriv, x, y, n);
        std::vector<double> res2 = Adams::Result(deriv, x, y, n);

        double h = 1. / abs(n);
        double err1 = 0, err2 = 0;
        for (size_t i = 0; i < res1.size(); i++)
        {
            double t1 = abs(res1[i] - answer(x)), t2 = abs(res2[i] - answer(x));
            
            if (t1 > err1)
                err1 = t1;
            if (t2 > err2)
                err2 = t2;

            x += h;
        }

        output << n << '\t' << std::setw(10) << h << '\t' << err1 << '\t' << err2 << std::endl;
    }
}