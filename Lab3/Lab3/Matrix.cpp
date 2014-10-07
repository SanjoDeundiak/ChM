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

Matrix::Matrix(size_t n_) :
	n(n_)
{	
	for (size_t i = 1; i <= n; i++)
	{
		vector<Item> T;
		Item t;
		cout << "Enter raw #" << i << ": ";
		for (size_t j = 1; j <= n; j++)
		{
			cin >> t;
			T.push_back(t);
		}
		A.push_back(T);
	}
	cout << endl;
}

void Matrix::LU(Matrix &L, Matrix &U) const
{
	U = *this;		

	for (size_t i = 0; i < n; i++)
	for (size_t j = i; j < n; j++)
	{
		if (abs(U.A[i][i]) < eps)
			throw std::logic_error("LU decompsition doesn't exist");
		L.A[j][i] = U.A[j][i] / U.A[i][i];		
	}	

	for (size_t k = 1; k < n; k++)
	{
		for (size_t i = k - 1; i < n; i++)
		for (size_t j = i; j < n; j++)
			L.A[j][i] = U.A[j][i] / U.A[i][i];

		for (size_t i = k; i < n; i++)
		for (size_t j = k - 1; j < n; j++)
			U.A[i][j] = U.A[i][j] - L.A[i][k - 1] * U.A[k - 1][j];
	}
}

void Matrix::Show() const
{
	std::streamsize pr = cout.precision(5);
	cout << std::left;
	cout << "MATRIX " << n << 'x' << n << ": " << endl;
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			cout << std::setw(10);
			cout << A[i][j];
			cout << ' ';
		}			
		cout << endl;
	}
	cout << endl;
	cout.precision(pr);
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

vector<Item> Matrix::Solve(const vector<Item> &b) const
{
	Matrix L(n, 0), U(n, 0);
	LU(L, U);	

	//(L*U).Show();

	return U.SolveU(L.SolveL(b));
}

Item Matrix::Determinant() const
{
	Matrix L(n, 0), U(n, 0);

	Item det = 1;

	LU(L, U);
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

vector<Item> operator+(const vector<Item> &v1, const vector<Item> &v2)
{
	size_t l = std::min(v1.size(), v2.size());
	vector<Item> r;
	r.resize(l);
	for (size_t i = 0; i < l; i++)
		r[i] = v1[i] + v2[i];

	return r;
}

void ShowV(const vector<double> &v)
{
	cout << "{";
	for (size_t i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << "}" << endl;
}

vector<Item> Matrix::Solve2(const vector<double> &b) const
{
	//if (!DiagDom())
		//throw std::exception("Matrix is not diagonally dominant");

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
		cout << "Iteration " << it << endl << "x=";
		ShowV(x);
		cout << "eps=";
		ShowV(e);
		cout << "|eps|=" << ep << endl << endl;
		x = L.SolveL(U*x + b);
		e = b - (*this)*x;
		ep = sqrt(e * e);
		++it;
		if (it > MAXIT)
			throw std::exception("Too many iterations");
	} 

	cout << "Result in " << it << " iterations:" << endl << "x=";
	ShowV(x);
	cout << "eps=";
	ShowV(e);
	cout << "|eps|=" << ep << endl;
	
	return x;
}