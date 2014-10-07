#include "Polynomial.h"
#include "Vector.h"
#include <iostream>

Polynomial::Polynomial(const vector<Item> &A) : vec(A)
{
}

Item Polynomial::operator()(Item x) const
{
	Item r = 0, p = 1;
	for (auto i = vec.begin(); i != vec.end(); i++)
	{
		r += *i * p;
		p *= x;
	}

	return r;
}

Polynomial Polynomial::operator+(const Polynomial &other) const
{
	return Polynomial(vec + other.vec);
}

Polynomial Polynomial::operator*(double n) const
{
	return Polynomial(vec * n);
}

std::ostream &operator<<(std::ostream &os, const Polynomial &p)
{	
	os.setf(std::ios::fixed, std::ios::floatfield);
	std::streamsize pr = os.precision(10);
	for (int i = p.vec.size() - 1; i >= 0; i--)
	if (abs(p.vec[i]) > eps)
	{
		os << std::showpos << p.vec[i] << std::noshowpos;
		if (i > 0)
			os << "x";
		if (i > 1)
			os << "^" << i;
		os << ' ';
	}
	os.precision(pr);
	
	return os;
}

Polynomial Polynomial::operator*(const Polynomial &other) const
{
	vector<Item> res(vec.size() + other.vec.size(), 0);
	for (size_t i = 0; i < vec.size(); i++)
		for (size_t j = 0; j < other.vec.size(); j++)
			res[i + j] += vec[i] * other.vec[j];

	if (!res.empty())
	{
		auto i = res.end() - 1;

		while (!res.empty() && abs(*i) < eps)
			i = res.erase(i) - 1;
	}
	
	return Polynomial(res);	
}

Polynomial Polynomial::Derivative() const
{
	vector<Item> der(1, 0);
	for (size_t i = 1; i < vec.size(); i++)
	{
		der[i - 1] = i * vec[i];
	}

	return Polynomial(der);
}