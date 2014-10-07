#pragma once

#include <fstream>

class Polynomial
{
	double m_a5, m_a4, m_a3, m_a2, m_a1, m_a0;
	double Value(double x) const;
	Polynomial Derivative() const;
	Polynomial();
	mutable std::ofstream debug;
	
	public:
		Polynomial(double a5, double a4, double a3, double a2, double a1, double a0, const std::string &debugfile = "debug.txt");
		~Polynomial();
		double GetRootBin(double left, double right, double acc) const;
		double GetRootChord(double left, double right, double acc) const;
		double GetRootNewt(double left, double right, double acc) const;
};