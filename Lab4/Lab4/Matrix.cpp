#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <algorithm>

const double eps = 1E-5;
const long MAXIT = 1000;

using std::cin;
//using std::cout;
using std::endl;
using std::ifstream;

extern std::ofstream cout;

Matrix::Matrix(const vector<vector<Item>> &v) :
	n(v.size())
{
	for (size_t i = 0; i < n; i++)
	{
		vector<Item> T;		
		for (size_t j = 0; j < n; j++)		
			T.push_back(v[i][j]);		
		A.push_back(T);
	}	
}

Matrix::Matrix(size_t n_, Item fill) : n(n_)
{
	A.resize(n);
	for (size_t i = 0; i < n; i++)
		A[i].push_back(fill);
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

Matrix Matrix::Inverse2() const
{
	if (n != 2)
		throw std::logic_error("n!=2");

	Matrix m_inv(n, 0);

	Item det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
	vector<Item> b, r;
	b.resize(n);
	r.resize(2);

	r[0] = A[1][1] / det;
	r[1] = -A[1][0] / det;

	m_inv.A[0] = r;

	r[0] = -A[0][1] / det;
	r[1] = A[0][0] / det;

	m_inv.A[1] = r;

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

Item operator*(const vector<Item> &v1, const vector<Item> &v2)
{
	size_t l = std::min(v1.size(), v2.size());
	Item r = 0;
	for (size_t i = 0; i < l; i++)
		r += v1[i] * v2[i];

	return r;
}