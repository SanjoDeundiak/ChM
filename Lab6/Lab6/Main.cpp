#include "Matrix.h"
#include "Vector.h"
#include <fstream>
#include "Logger.h"

logger::log GlobalLog("log.txt");

int main()
{
    std::ifstream input("input.txt");
    Matrix m(4, input);

    /*{
        std::vector< std::vector<double>> eigenVectors;
        std::vector<double> eigenValues;
        m.eigenValuesJacobi(eigenValues, eigenVectors);
    }*/
    {
        double eigenValue = 0;
        std::vector<double> eigenVector;
        m.eigenValuesScalar(eigenValue, eigenVector);

        Matrix b(m);
        for (int i = 0; i < 4; i++)
            b[i][i] -= eigenValue;
        double newEigenValue = 0;
        eigenVector.clear();
        b.eigenValuesScalar(newEigenValue, eigenVector);
        newEigenValue += eigenValue;
        GlobalLog << "FINAL RESULT = " << newEigenValue << std::endl;
    }
    /*{
        std::vector< std::vector<double>> eigenVectors;
        std::vector<double> eigenValues;
        m.eigenValuesDanilevski(eigenValues, eigenVectors);
    }*/
}