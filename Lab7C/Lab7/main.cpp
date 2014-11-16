// y' = (2x + y - 3) / (x - 1)
// y() = 

#include "RKM.h"
#include <iostream>
#include <functional>
#include <vector>

int main()
{
    std::function<double(double, double)> deriv = [](double x, double y)->double{ return (2 * x + y - 3) / (x - 1); };

    std::function<double(double)> result = RKM::Solution(deriv, 2, 1);
    std::vector<double> res = RKM::Result(deriv, 2, 1);

    // std::cout << result(3) << std::endl;

    double x = 2, h = 0.1;
    for (double y : res)
    {
        std::cout << x << ' ' << y << std::endl;
        x += h;
    }

    std::cin.get();
}