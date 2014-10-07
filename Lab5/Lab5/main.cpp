#include "Interpolation.h"
#include <iostream>
#include <fstream>

using std::endl;

double F(double x)
{
	return pow(log(x + cos(x)), 2);
}
//
//int main()
//{
//	Interpolation inter;
//
//	Polynomial p({ 0 });
//	Spline s(0);
//
//	try
//	{
//		p = inter.Lagrange();
//		s = inter.Spline2();
//		std::cout << "Lagrange polynomial:\n" << inter.Lagrange() << endl << endl
//				<< "Newton forawrd polynomial\n" << inter.Newton1() << endl << endl
//				<< "Newton backward polynomial\n" << inter.Newton2() << endl << endl << endl
//				<< "Quadratic spline:\n" << inter.Spline2() << endl;
//	}
//	catch (std::exception e)
//	{
//		std::cout << e.what();
//	}
//
//	std::ofstream debug("debug.txt");
//	{
//		debug << "Precision analysis:\nNewton polynomial\n";
//		Item x = inter.Left(), dx = (inter.Right() - inter.Left()) / (inter.n() - 1) / 5;
//		while (x - inter.Right() < eps)
//		{
//			debug.precision(5);
//			debug.setf(std::ios::fixed, std::ios::floatfield);
//			debug << "x=" << x << " y=" << F(x) << " Interpol=" << p(x) << " delta=" << abs(F(x) - p(x)) << endl;
//			x += dx;
//		}
//	}
//
//	{
//		debug << "\nQuadratic spline polynomial\n";
//		Item x = inter.Left(), dx = (inter.Right() - inter.Left()) / (inter.n() - 1) / 5;
//		while (x - inter.Right() < eps)
//		{
//			debug.precision(5);
//			debug.setf(std::ios::fixed, std::ios::floatfield);
//			debug << "x=" << x << " y=" << F(x) << " Interpol=" << s(x) << " delta=" << abs(F(x) - s(x)) << endl;
//			x += dx;
//		}
//	}
//	std::cin.get();
//}