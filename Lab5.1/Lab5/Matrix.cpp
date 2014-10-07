#include "Matrix.h"
#include <iomanip>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <algorithm>
#include "Vector.h"

extern const double eps;
const long MAXIT = 1000;

using std::cin;
//using std::cout;
using std::endl;
using std::ifstream;

Matrix::Matrix(size_t n_, Item fill) :
	n(n_)
{
	for (size_t i = 1; i <= n; i++)
	{
		vector<Item> T;		
		for (size_t j = 1; j <= n; j++)		
			T.push_back(fill);		
		A.push_back(T);
	}	
}

Matrix::Matrix(const vector<vector<Item>> &A_) : A(A_), n(A.size())
{
	size_t n = A.size();
	for (size_t i = 0; i < n; i++)
	if (A[i].size() != n)
		throw std::invalid_argument("Not square matrix");
}

Matrix::Matrix(size_t n_, ifstream &file) :
	n(n_)
{	
	for (size_t i = 1; i <= n; i++)
	{
		vector<Item> T;
		Item t;		
		for (size_t j = 1; j <= n; j++)
		{
			file >> t;
			T.push_back(t);
		}
		A.push_back(T);
	}
}

void Matrix::LU(Matrix &L, Matrix &U, vector<Item> &v) const
{	
	U = *this;
	for (size_t i = 0; i < n; i++) // for all diagonal elements
	{
		{ // find max element in raw and place it in diagonal position
			Item max = abs(U.A[i][i]);
			size_t maxi = i;
			for (size_t j = i; j < n; j++)
			if (abs(U.A[j][i]) > max)
			{
				max = abs(U.A[j][i]);
				maxi = j;
			}
			if (max < eps)
				throw std::runtime_error("LU decomposition doesn't exists");
			std::swap(U.A[i], U.A[maxi]);
			std::swap(L.A[i], L.A[maxi]);
			std::swap(v[i], v[maxi]);
		}
		for (size_t j = i; j < n; j++)
		{
			L.A[j][i] = U.A[j][i] / U.A[i][i];
			if (i !=j)
				U.A[j] = U.A[j] + (-L.A[j][i]) * U.A[i];
		}		
	}	
}

Matrix Matrix::operator*(const Matrix &other) const
{
	Matrix P(n, 0);

	for (size_t row = 0; row < n; row++) 	
		for (size_t col = 0; col < n; col++) 		
			for (size_t inner = 0; inner < n; inner++) 			
				P.A[row][col] += A[row][inner] * other.A[inner][col];

	return P;
}

vector<Item> Matrix::SolveL(const vector<Item> &b) const
{	
	vector<Item> r;
	r.reserve(n);
	Item t;
	for (size_t i = 0; i < n; i++)
	{
		t = 0;
		for (size_t j = 0; j < i; j++)
			t += A[i][j] * r[j];
		r.push_back((b[i] - t) / A[i][i]);
	}

	return r;
}

vector<Item> Matrix::SolveU(const vector<Item> &b) const
{
	vector<Item> r;
	r.resize(n);
	Item t;
	for (size_t i = n; i > 0; i--)
	{
		t = 0;
		for (size_t j = n-1; j > i - 1; j--)
			t += A[i - 1][j] * r[j];
		r[i - 1] = (b[i - 1] - t) / A[i - 1][i - 1];
	}

	return r;
}

vector<Item> Matrix::Solve(vector<Item> b) const
{
	Matrix L(n, 0), U(n, 0);
	LU(L, U, b);		

	return U.SolveU(L.SolveL(b));
}

Item Matrix::Determinant() const
{
	Matrix L(n, 0), U(n, 0);

	Item det = 1;

	LU(L, U, vector<Item>(n, 0));
	for (size_t i = 0; i < n; i++)
		det *= U.A[i][i];

	return det;
}

Matrix Matrix::Inverse()
{
	Matrix m_inv(n, 0);

	vector<Item> b, r;
	b.resize(n);

	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
			b[j] = (i == j) ? 1 : 0;
		r = Solve(b);

		for (size_t j = 0; j < n; j++)
			(m_inv.A)[j][i] = r[j];
	}

	return m_inv;
}

vector<Item> Matrix::operator*(const vector<Item> &vec) const
{
	vector<Item> res;
	for (size_t i = 0; i < n; i++)
	{
		Item r = 0;
		for (size_t j = 0; j < n; j++)
			r += A[i][j] * vec[j];
		res.push_back(r);
	}

	return res;
}

bool Matrix::DiagDom() const
{
	Item s;
	for (size_t i = 0; i < n; i++)
	{
		s = 0;
		for (size_t j = 0; j < n; j++)
			s += (i != j) ? abs(A[i][j]) : 0.0;
		if (abs(A[i][i]) < s)
			return false;
	}

	return true;
}

Item operator*(const vector<Item> &v1, const vector<Item> &v2)
{
	size_t l = std::min(v1.size(), v2.size());
	Item r = 0;
	for (size_t i = 0; i < l; i++)
		r += v1[i] * v2[i];

	return r;
}

vector<Item> operator-(const vector<Item> &v1, const vector<Item> &v2)
{
	size_t l = std::min(v1.size(), v2.size());
	vector<Item> r;
	r.resize(l);
	for (size_t i = 0; i < l; i++)
		r[i] = v1[i] - v2[i];

	return r;
}

vector<Item> Matrix::Solve2(const vector<double> &b) const
{		
	Matrix L = *this, U = *this;
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = i + 1; j < n; j++)
			L.A[i][j] = 0;
	}
	for (size_t i = 0; i < n; i++)
	for (size_t j = 0; j < n; j++)
		U.A[i][j] = (j <= i) ? 0 : -A[i][j];

	vector<Item> x = b, e = b - (*this)*x;;
	Item ep = sqrt(e*e);
	long it = 1;
	
	while (ep > eps)
	{				
		x = L.SolveL(U*x + b);
		e = b - (*this)*x;
		ep = sqrt(e * e);
		++it;
		if (it > MAXIT)
			throw std::exception("Too many iterations");
	} 
	
	return x;
}

std::ostream &operator<<(std::ostream &os, const Matrix &m)
{
	os << std::left;
	for (size_t i = 0; i < m.A.size(); i++)
	{
		for (size_t j = 0; j < m.A[i].size(); j++)
			os << std::setw(9) << m.A[i][j] << " ";

		os << std::endl;
	}

	return os;
}