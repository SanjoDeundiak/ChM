#include "Polynomial.h"
#include "Vector.h"
#include <algorithm>

vector<Item> operator+(const vector<Item> &v1, const vector<Item> &v2)
{	
	size_t lmin = v1.size(), lmax = v2.size();

	if (lmin > lmax)
		std::swap(lmin, lmax);
	
	vector<Item> res(lmax);

	for (size_t i = 0; i < lmin; i++)
		res[i] = v1[i] + v2[i];

	for (size_t i = lmin; i < lmax; i++)
		res[i] = (lmax == v1.size()) ? v1[i] : v2[i];

	return res;
}

vector<Item> operator*(const vector<Item> v, double n)
{
	vector<Item> r(v.size());
	for (size_t i = 0; i < r.size(); i++)
		r[i] = v[i] * n;

	return r;
}

vector<Item> operator*(double n, const vector<Item> v)
{
	return v * n;
}

std::ostream &operator<<(std::ostream &os, const vector<Item> &v)
{
	for (size_t i = 0; i < v.size(); i++)
		os << v[i] << " ";	

	return os;
}