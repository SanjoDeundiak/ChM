#include "Spline.h"
#include <algorithm>
#include <iostream>

Spline::Spline(Item x)
{
	X.push_back(x);
}

void Spline::Add(const Polynomial &p, Item x)
{
	P.push_back(p);
	X.push_back(x);
}

std::ostream &operator<<(std::ostream &os, const Spline &s)
{
	for (size_t i = 0; i < s.P.size(); i++)
	{
		os << s.P[i] << " on [" << s.X[i] << "," << s.X[i + 1] << "]\n";
	}	
	
	return os;
}

Item Spline::operator()(Item x) const
{	
	size_t i;
	for (i = 0; i < P.size() - 1 && x > X[i + 1]; i++);

	return P[i](x);
}
