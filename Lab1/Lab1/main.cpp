#include "Polynomial.h"
#include <iostream>

const double eps = 1E-05;

using std::endl;
using std::cout;

int main()
{
	Polynomial pol(0, 1, -2, -9, -3, -3);
	//Polynomial pol(16, -3, 1, 2, -4, 2);
	
	double a = -3.5, b = -1.8;

	cout.setf(std::ios::scientific);
	try
	{
		cout << pol.GetRootBin(a, b, eps) << endl;
	}
	catch (std::exception e)
	{
		cout << e.what() << endl;
	}
	
	try
	{
		cout << pol.GetRootChord(a, b, eps) << endl;
	}
	catch (std::exception e)
	{
		cout << e.what() << endl;
	}

	try
	{
		cout << pol.GetRootNewt(a, b, eps) << endl;
	}	
	catch (std::exception e)
	{
		cout << e.what() << endl;
	}

	cout.unsetf(std::ios::scientific);


	std::cin.get();
	std::cin.get();

	return 0;
}