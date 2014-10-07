#include "Interpolation.h"
#include "Vector.h"
#include "Polynomial.h"
#include "Spline.h"
#include "Matrix.h"
#include <fstream>
#include <iostream>
#include <map>
//#include <cmath>

const double eps = 1E-6;

const std::string file = "input.txt";

Interpolation::Interpolation()
{
	std::ifstream cin(file);
	int n;
	if (!(cin >> n))
		throw std::runtime_error("Error reading from the file");

	for (int i = 0; i < n && cin.good(); i++)
	{
		Item x, y;
		cin >> x >> y;
		X_.push_back(x);
		Y_.push_back(y);
	}	
}

Interpolation::Interpolation(Item a, Item b, int n)
{
	Item x = a, dx = (b - a) / n;

	for (int i = 0; i < n; i++)
	{
		X_.push_back(x);
		Y_.push_back(OriginalF(x));
		x += dx;
	}
}

Item Interpolation::OriginalF(double x) const
{
	return pow(log(x + cos(x)), 2);
}

Polynomial Interpolation::Lagrange() const
{
	return Lagrange_(X_, Y_);
}


Polynomial Interpolation::Lagrange_(const vector<Item> &X, const vector<Item> &Y) const
{
	size_t n = X.size();
	Polynomial L({ 0 });
	for (size_t i = 0; i < n; i++)
	{
		Polynomial l({ 1 });
		for (size_t j = 0; j < n; j++)		
		if (i != j)
				l = l * Polynomial({ -X[j] / (X[i] - X[j]), 1 / (X[i] - X[j])});		

		L = L + Y[i] * l;
	}

	return L;
}

Polynomial Interpolation::Newton1() const
{
	return Newton1_(X_, Y_);
}

Polynomial Interpolation::Newton2() const
{
	vector<Item> X1(X_.size()), Y1(X_.size());

	for (size_t i = 0; i < X1.size(); i++)
	{
		X1[i] = X_[X_.size() - i - 1];
		Y1[i] = Y_[X_.size() - i - 1];
	}

	return Newton1_(X1, Y1);
}


Polynomial Interpolation::Newton1_(const vector<Item> &X, const vector<Item> &Y) const
{
	size_t n = X.size(), i = 0;
	vector<vector<Item>> Z(n, vector<Item>(n));		
	for (size_t i = 0; i < n; i++)
		Z[i][0] = Y[i];			
		
	for (size_t i = 1; i < n; i++)
	{			
		size_t it1 = i;			
		for (size_t j = 0; j < n - i; j++, it1++)
			Z[j][i] = (Z[j][i - 1] - Z[j + 1][i - 1]) / (X[j] - X[it1]);
	}	

	Polynomial L({ Z[0][0] }), l({ 1 });	
	for (size_t i = 1; i < n; i++)
	{
		l = l * Polynomial({-X[i - 1], 1});
		L = L + Z[0][i] * l;
	}

	return L;
}

Spline Interpolation::Spline2() const
{
	const int c = 3;
	size_t n = X_.size(), k = c * (n - 1);
	vector<vector<Item>> A(k, vector<Item>(k, 0));
	vector<Item> b(k, 0);
	Spline S(X_[0]);
	for (size_t i = 0; i < n - 1; i++) // values
	{
		for (size_t j = 0; j < c; j++)
		{
			A[i * 2][c * i + j] = pow(X_[i], j);
			A[i * 2 + 1][c * i + j] = pow(X_[i + 1], j);
		}
		b[i * 2] = Y_[i];
		b[i * 2 + 1] = Y_[i + 1];
	}

	for (size_t i = 1; i < n - 1; i++) // derivatives
	{
		for (size_t j = 1; j < c; j++)
		{
			double a = pow(X_[i], j - 1), b = j * a;
			A[i + 2 * n - 3][j + c*(i - 1)] = j * pow(X_[i], j - 1);
		}
		for (size_t j = 1; j < c; j++)		
			A[i + 2*n - 3][j + c * i] = -int(j) * pow(X_[i], j - 1);		
		b[i + 2*n - 3] = 0;
	}

	for (size_t j = 1; j < c; j++) // edge condition
		A[3 * n - 4][j] = j * pow(X_[0], j - 1);
	b[3 * n - 4] = 0;
	
	Matrix m(A);	
	vector<Item> r = m.Solve(b), r1;		
	r1.resize(c);
	for (size_t i = 0; i < n - 1; i++)
	{
		for (size_t j = 0; j < c; j++)
			r1[j] = r[i * c + j];
		S.Add(Polynomial(r1), X_[i + 1]);
	}	
	return S;
}
