#include "Polynomial.h"
#include "cmath"
#include <string>

using std::endl;

template <typename T>
inline bool sgn(T val) 
{
	return (val >= T(0));
}

Polynomial::Polynomial(double a5, double a4, double a3, double a2, double a1, double a0, const std::string &debugfile) :
	m_a5(a5), m_a4(a4),
	m_a3(a3), m_a2(a2),
	m_a1(a1), m_a0(a0)
{
	debug.open(debugfile, std::ios_base::out);
	debug.precision(7);

}

Polynomial::~Polynomial()
{
	debug.close();
}

inline double Polynomial::Value(double x) const
{
	return m_a5*pow(x, 5) + m_a4*pow(x, 4) + m_a3*pow(x, 3) + m_a2*pow(x, 2) + m_a1*x + m_a0;
}

double Polynomial::GetRootBin(double left, double right, double acc) const
{
	debug << "Bisection method" << endl;

	long it = 0;
	double mid, midvalue;

	if (sgn(Value(left)) == sgn(Value(right)))
		throw std::runtime_error("No root or there are more then 1");

	do
	{
		++it;
		mid = (left + right) / 2;
		midvalue = Value(mid);
		if (midvalue * Value(left) > 0)		
			left = mid;					
		else		
			right = mid;
		debug << "Iteration " << it << " [" << left << ", " << right << "]"
			<< " f1=" << Value(left) << " f2=" << Value(right) << endl;
	} while (right - left > acc || abs(midvalue) > acc);

	debug << "Result: x = " << mid << " in " << it << " iterations" << endl << endl;

	return mid;
}

double Polynomial::GetRootChord(double left, double right, double acc) const
{	
	debug << "Chord method" << endl;
	long it = 0;
	double c0, c = 0, cvalue;	

	if (sgn(Value(left)) == sgn(Value(right)))
		throw std::runtime_error("No root or there are more then 1");

	do
	{
		++it;
		c0 = c;
		c = (left * Value(right) - right * Value(left)) / (Value(right) - Value(left));
		cvalue = Value(c);

		if (sgn(Value(left)) != sgn(cvalue))
			right = c;
		else
			left = c;
		debug << "Iteration " << it << " [" << left << ", " << right << "]" 
			<< " f1=" << Value(left) << " f2=" << Value(right) << endl;
	} while (abs(c - c0) > acc || abs(cvalue) > acc);

	debug << "Result: x = " << c << " in " << it << " iterations" << endl << endl;

	return c;
}

double Polynomial::GetRootNewt(double left, double right, double acc) const
{
	debug << "Newton's method" << endl;
	long it = 1;
	Polynomial derivative = Derivative();		
	double c0, c, cvalue;	

	if (sgn(Value(left)) == sgn(Value(right)))
		throw std::runtime_error("No root or there are more then 1");

	c = left - Value(left) / derivative.Value(left);	
	c0 = right - Value(right) / derivative.Value(right);
	
	if (abs(Value(c0)) < abs(Value(c)))
	{
		c = c0;
	}

	cvalue = Value(c);

	debug << "Iteration " << 1 << " x=" << c << " f(x)=" << cvalue << endl;

	while (abs(c - c0) > acc || abs(cvalue) > acc)
	{			
		++it;
		c0 = c;
		c = c0 - Value(c) / derivative.Value(c);
		cvalue = Value(c);
		debug << "Iteration " << it << " x=" << c << " f(x)=" << cvalue << endl;
	}

	debug << "Result: x = " << c << " in " << it << " iterations" << endl << endl;
	return c;
}

Polynomial Polynomial::Derivative() const
{
	return Polynomial(0, 5 * m_a5, 4 * m_a4, 3 * m_a3, 2 * m_a2, m_a1, "debug1.txt");
}