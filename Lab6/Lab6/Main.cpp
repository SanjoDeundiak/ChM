#include "Matrix.h"
#include "Vector.h"
#include <fstream>
#include "Logger.h"

logger::log GlobalLog("log.txt");

int main()
{
    std::ifstream input("input.txt");
    Matrix m(4, input);
    std::vector< std::vector<double>> eigenVectors;
    std::vector<double> eigenValues;

    double eigenValue = 0;
    std::vector<double> eigenVector;
   
    //m.eigenValuesJacobi(eigenValues, eigenVectors);
    //m.eigenValuesScalar(eigenValue, eigenVector);
    m.eigenValuesDanilevski(eigenValues, eigenVectors);
}