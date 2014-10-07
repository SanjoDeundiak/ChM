#include "Vector.h"
#include <iostream>
#include <algorithm>

using std::vector;

vector<double> operator+(const vector<double> &v1, const vector<double> &v2)
{	
	size_t lmin = v1.size(), lmax = v2.size();

	if (lmin > lmax)
		std::swap(lmin, lmax);
	
	vector<double> res(lmax);

	for (size_t i = 0; i < lmin; i++)
		res[i] = v1[i] + v2[i];

	for (size_t i = lmin; i < lmax; i++)
		res[i] = (lmax == v1.size()) ? v1[i] : v2[i];

	return res;
}

vector<double> operator*(const vector<double> v, double n)
{
	vector<double> r(v.size());
	for (size_t i = 0; i < r.size(); i++)
		r[i] = v[i] * n;

	return r;
}

vector<double> operator*(double n, const vector<double> v)
{
	return v * n;
}

std::ostream &operator<<(std::ostream &os, const vector<double> &v)
{
	for (size_t i = 0; i < v.size(); i++)
		os << v[i] << " ";	

	return os;
}