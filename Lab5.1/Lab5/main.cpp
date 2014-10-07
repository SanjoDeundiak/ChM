#include "Interpolation.h"
#include <iostream>
#include <fstream>

using std::endl;

double F(double x)
{
	return pow(log(x + cos(x)), 2);
}

int main()
{
	Interpolation inter;

	Polynomial p1({ 0 }), p2({ 0 }), p3({ 0 });
	Spline s2(0), s3(0);

	try
	{
		p1 = inter.Lagrange(); p2 = inter.Newton1(); p3 = inter.Newton2();
		s2 = inter.Spline2(); s3 = inter.Spline3();
		std::cout << "Lagrange polynomial:\n" << p1 << endl << endl
				<< "Newton forawrd polynomial\n" << p2 << endl << endl
				<< "Newton backward polynomial\n" << p3 << endl << endl << endl
				<< "Quadratic spline:\n" << s2 << endl << endl
		        << "Qubic spline:\n" << s3 << endl;
	}
	catch (std::exception e)
	{
		std::cout << e.what();
	}

	std::ofstream debug("debug.txt");
	{
		debug << "Precision analysis:\nNewton polynomial\n";
		Item x = inter.Left(), dx = (inter.Right() - inter.Left()) / (inter.n() - 1) / 5;
		while (x - inter.Right() < eps)
		{
			debug.precision(5);
			debug.setf(std::ios::fixed, std::ios::floatfield);
			debug << "x=" << x << " y=" << F(x) << " Interpol=" << p1(x) << " delta=" << abs(F(x) - p1(x)) << endl;		
			x += dx;
		}
	}

	{
		debug << "\nQuadratic spline polynomial\n";
		Item x = inter.Left(), dx = (inter.Right() - inter.Left()) / (inter.n() - 1) / 5;
		while (x - inter.Right() < eps)
		{
			debug.precision(5);
			debug.setf(std::ios::fixed, std::ios::floatfield);
			debug << "x=" << x << " y=" << F(x) << " Interpol=" << s2(x) << " delta=" << abs(F(x) - s2(x)) << endl;
			x += dx;
		}
	}

	{
		debug << "\nQubic spline polynomial\n";
		Item x = inter.Left(), dx = (inter.Right() - inter.Left()) / (inter.n() - 1) / 5;
		while (x - inter.Right() < eps)
		{
			debug.precision(5);
			debug.setf(std::ios::fixed, std::ios::floatfield);
			//debug << "x=" << x << " y=" << F(x) << " Interpol=" << s3(x) << " delta=" << abs(F(x) - s3(x)) << endl;
			debug << abs(F(x) - s3(x)) << endl;
			x += dx;
		}
	}

	{
		Polynomial m({4.18});
		m = m * Polynomial({0, 1});
		m = m * Polynomial({ -1, 1 });
		m = m * Polynomial({ -2, 1 });
		m = m * Polynomial({ -3, 1 });
		m = m * Polynomial({ -4, 1 });
		debug << "\Majorant polynomial\n";
		Item x = inter.Left(), dx = (inter.Right() - inter.Left()) / (inter.n() - 1) / 5;
		while (x - inter.Right() < eps)
		{
			debug.precision(5);
			debug.setf(std::ios::fixed, std::ios::floatfield);
			debug << "x=" << x << " y=" << F(x) << " Major=" << abs(m(x)) << endl;
			x += dx;
		}
	}
	std::cin.get();
}