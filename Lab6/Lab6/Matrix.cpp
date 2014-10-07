#include "Matrix.h"
#include <iomanip>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <algorithm>
#include "Vector.h"
#include "Logger.h"
#include "Polynomial.h"

extern logger::log GlobalLog;

#define M_PI 3.1415926535897932384626433832795

static const double eps = 1E-05;
const long MAXIT = 1000;

using std::cin;
//using std::cout;
using std::endl;
using std::ifstream;

using std::vector;
using std::string;
using std::ifstream;

Matrix::Matrix(size_t n_, double fill) :
	n(n_)
{
	for (size_t i = 1; i <= n; i++)
	{
		vector<double> T;		
		for (size_t j = 1; j <= n; j++)		
			T.push_back(fill);		
		A.push_back(T);
	}	
}

Matrix::Matrix(size_t dimension, bool identity /*=false*/) :
    n(dimension)
{
    A.resize(n);
    for (size_t i = 0; i < n; i++)
    {
        A[i].resize(n);
        if (identity)
        for (size_t j = 0; j < n; j++)
            A[i][j] = (i == j) ? 1 : 0;
    }
}

Matrix::Matrix(const vector<vector<double>> &A_) : A(A_), n(A.size())
{
	size_t n = A.size();
	for (size_t i = 0; i < n; i++)
	if (A[i].size() != n)
		throw std::invalid_argument("Not square matrix");
}

Matrix::Matrix(size_t dimension, std::ifstream &file) :
    n(dimension)
{	
	for (size_t i = 1; i <= n; i++)
	{
		vector<double> T;
		double t;		
		for (size_t j = 1; j <= n; j++)
		{
			file >> t;
			T.push_back(t);
		}
		A.push_back(T);
	}
}

std::vector<double> Matrix::GetColumn(size_t i) const
{
    std::vector<double> result(n);
    for (size_t j = 0; j < n; j++)
    {
        result[j] = A[j][i];
    }
    return result;
}

void Matrix::LU(Matrix &L, Matrix &U, vector<double> &v) const
{	
	U = *this;
	for (size_t i = 0; i < n; i++) // for all diagonal elements
	{
		{ // find max element in raw and place it in diagonal position
			double max = abs(U.A[i][i]);
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
	Matrix P(n, 0.0);

	for (size_t row = 0; row < n; row++) 	
		for (size_t col = 0; col < n; col++) 		
			for (size_t inner = 0; inner < n; inner++) 			
				P.A[row][col] += A[row][inner] * other.A[inner][col];

	return P;
}

vector<double> Matrix::SolveL(const vector<double> &b) const
{	
	vector<double> r;
	r.reserve(n);
	double t;
	for (size_t i = 0; i < n; i++)
	{
		t = 0;
		for (size_t j = 0; j < i; j++)
			t += A[i][j] * r[j];
		r.push_back((b[i] - t) / A[i][i]);
	}

	return r;
}

vector<double> Matrix::SolveU(const vector<double> &b) const
{
	vector<double> r;
	r.resize(n);
	double t;
	for (size_t i = n; i > 0; i--)
	{
		t = 0;
		for (size_t j = n-1; j > i - 1; j--)
			t += A[i - 1][j] * r[j];
		r[i - 1] = (b[i - 1] - t) / A[i - 1][i - 1];
	}

	return r;
}

vector<double> Matrix::Solve(vector<double> b) const
{
	Matrix L(n, 0.0), U(n, 0.0);
	LU(L, U, b);		

	return U.SolveU(L.SolveL(b));
}

double Matrix::Determinant() const
{
	Matrix L(n, 0.0), U(n, 0.0);

	double det = 1;

	LU(L, U, vector<double>(n, 0));
	for (size_t i = 0; i < n; i++)
		det *= U.A[i][i];

	return det;
}

Matrix Matrix::Inverse()
{
	Matrix m_inv(n, 0.0);

	vector<double> b, r;
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

vector<double> Matrix::operator*(const vector<double> &vec) const
{
	vector<double> res;
	for (size_t i = 0; i < n; i++)
	{
		double r = 0;
		for (size_t j = 0; j < n; j++)
			r += A[i][j] * vec[j];
		res.push_back(r);
	}

	return res;
}

bool Matrix::DiagDom() const
{
	double s;
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

double operator*(const vector<double> &v1, const vector<double> &v2)
{
	size_t l = std::min(v1.size(), v2.size());
	double r = 0;
	for (size_t i = 0; i < l; i++)
		r += v1[i] * v2[i];

	return r;
}

vector<double> operator-(const vector<double> &v1, const vector<double> &v2)
{
	size_t l = std::min(v1.size(), v2.size());
	vector<double> r;
	r.resize(l);
	for (size_t i = 0; i < l; i++)
		r[i] = v1[i] - v2[i];

	return r;
}

vector<double> Matrix::Solve2(const vector<double> &b) const
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

	vector<double> x = b, e = b - (*this)*x;;
	double ep = sqrt(e*e);
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
			os << std::setw(15) << m.A[i][j] << " ";

		os << std::endl;
	}

	return os;
}

void Matrix::findMax(size_t &i, size_t &j) const
{
    double max = 0;
    for (size_t _i = 0; _i < n; _i++)
    for (size_t _j = _i + 1; _j < n; _j++)
    if (abs(A[_i][_j]) > max)
    {
        i = _i;
        j = _j;
        max = abs(A[i][j]);
    }
}

Matrix Matrix::Trans() const
{
    Matrix temp(n);

    for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
        temp.A[i][j] = A[j][i];

    return temp;
}

double Matrix::sigma() const
{
    double sum = 0;
    for (size_t i = 0; i < n; i++)
        sum += pow(A[i][i], 2);

    return sum;
}

double Matrix::omega() const
{
    double sum = 0;
    for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < i; j++)
        sum += pow(A[i][j], 2);

    return sum;
}

size_t Matrix::eigenValuesJacobi(vector<double> &eigenValues, vector< vector<double> > &eigenVector, size_t itLimit /*= 1000*/) const
{
    GlobalLog << "####JACOBI METHOD####" << std::endl << std::endl;

    Matrix Ai(*this), B(n, true);

    size_t it;
    for (it = 0; it < itLimit; it++)
    {
        GlobalLog << "Iteration #" << it << std::endl;
        GlobalLog << "Matrix:" << std::endl << Ai << std::endl;

        size_t i, j;
        Ai.findMax(i, j);

        double fi = (abs(Ai[i][i] - Ai[j][j]) > eps) ? 0.5 * atan(2 * Ai[i][j] / (Ai[i][i] - Ai[j][j])) : M_PI / 4;

        GlobalLog << "Maximal element is: " << Ai[i][j] << " which is in position (" << i << "," << j << ")" << std::endl << std::endl;
        GlobalLog << "alpha=a[i][i]=" << Ai[i][i] << std::endl
            << "beta=a[j][j]=" << Ai[j][j] << std::endl
            << "gamma=a[i][j]=" << Ai[i][j] << std::endl << std::endl;
        GlobalLog << "fi=" << fi << std::endl
            << "cos(fi)=" << cos(fi) << std::endl
            << "sin(fi)=" << sin(fi) << std::endl << std::endl;
        GlobalLog << "sigma=" << Ai.sigma() << std::endl
            << "omega=" << Ai.omega() << std::endl
            << "sum=" << Ai.sigma() + 2 * Ai.omega() << std::endl << std::endl;

        if (pow(Ai.omega(), 0.5) < eps)
            break;

        Matrix U(n, 0.0);

        for (size_t _i = 0; _i < n; _i++)
            U[_i][_i] = (_i != i && _i != j) ? 1 : 0;

        U[i][i] = cos(fi);
        U[j][j] = cos(fi);
        U[i][j] = -sin(fi);
        U[j][i] = sin(fi);

        // iteration
        Ai = U.Trans() * Ai * U;
        B = B * U;
    }

    eigenValues.clear();
    for (size_t i = 0; i < n; i++)
    {
        eigenValues.push_back(Ai[i][i]);
        eigenVector.push_back(B.GetColumn(i));
    }

    GlobalLog << "---RESULTS---" << std::endl
        << "Finished in " << it << " iterations" << std::endl << std::endl;

    for (size_t i = 0; i < eigenValues.size(); i++)
    {
        GlobalLog << "Eigen value #" << i + 1 << ": Value=" << eigenValues[i] << "\tVector: " << eigenVector[i] << std::endl
            << "Error=" << *this * eigenVector[i] - eigenValues[i] * eigenVector[i] << std::endl << std::endl;
    }

    return it;
}

size_t Matrix::eigenValuesScalar(double eigenValue, std::vector<double> &eigenVector, size_t itLimit) const
{
    GlobalLog << "####SCALAR METHOD####" << std::endl << std::endl;

    vector<double> &y = eigenVector;
    y.clear();
    y = vector<double>(n, 1);
    eigenValue = 0;
    double old = 0;

    size_t it;
    for (it = 0; it < itLimit; it++)
    {
        GlobalLog << "Iteration #" << it << std::endl;
        GlobalLog << "lambda=" << eigenValue << " y=" << y << std::endl << std::endl;

        if (abs(old - eigenValue) < eps && it > 1)
            break;

        old = eigenValue;
        eigenValue = (((*this * y) * y) / (y * y));
        y = *this * y;
        y = y * pow(y * y, -0.5);
    }

    GlobalLog << "---RESULTS---" << std::endl
        << "Finished in " << it << " iterations" << std::endl << std::endl;

    GlobalLog << "Eigen value=" << eigenValue << "\tVector: " << eigenVector << std::endl
        << "Error=" << *this * eigenVector - eigenValue * eigenVector << std::endl << std::endl;

    return it;
}

void Matrix::CCInvers(Matrix& B, Matrix& C) const
{
    for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
    {
        if (j == n - 1)
            B[i][j] = A[i][j];
        else if (i == j + 1)
            B[i][j] = 1;
        else
            B[i][j] = 0;

        if (j == 0)
            C[i][j] = ((i < n - 1) ? -A[i + 1][n - 1] : 1) / A[0][n - 1];
        else if (i == j - 1)
            C[i][j] = 1;
        else
            C[i][j] = 0;
    }
}

void Matrix::eigenValuesDanilevski(std::vector<double> &eigenValuesvector, std::vector< std::vector<double> > &eigenVector) const
{
    GlobalLog << "####DANILEVSKI METHOD####" << std::endl << std::endl;

    Matrix B(*this);
    Matrix C1(n), C2(n);

    for (size_t i = 0; i < n - 1; i++)
    {
        B.CCInvers(C1, C2);
        GlobalLog << "ITERATION #" << i << std::endl;
        GlobalLog << "C\n" << C1 << std::endl << std::endl;
        GlobalLog << "C-1\n" << C2 << std::endl << std::endl;
        B = C2 * B * C1;
        GlobalLog << "B\n" << B << std::endl << std::endl;
    }

    std::vector<double> pol = reverse(B.GetColumn(n - 1));
    pol = pol * (-1);
    pol.push_back(1);

    Polynomial p(pol);

    GlobalLog << "RESULT: " << std::endl
        << p;
}

