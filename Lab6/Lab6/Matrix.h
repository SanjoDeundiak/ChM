#pragma once

#include <vector>
#include <string>
#include <fstream>

class Matrix
{
    std::vector< std::vector<double> > A;
	size_t n;

    std::vector<double> SolveL(const std::vector<double> &b) const;
    std::vector<double> SolveU(const std::vector<double> &b) const;

    void findMax(size_t &i, size_t &j) const;

    std::vector<double> &operator[](size_t i) { return A[i]; }

    void CCInvers(Matrix& B, Matrix& C) const;
public:
	Matrix(size_t dimension, double fill);
	Matrix(size_t dimension, bool identity = false);
    Matrix(size_t dimension, std::ifstream &file);
    Matrix(const std::vector< std::vector<double> > &A);

    std::vector<double> GetColumn(size_t i) const;

	double Determinant() const;
	bool DiagDom() const;

    Matrix Inverse();

    Matrix operator*(const Matrix &other) const;
    std::vector<double> operator*(const std::vector<double> &vec) const;

    const std::vector<double> &operator[](size_t i) const { return A[i]; }

    std::vector<double> Solve(std::vector<double> b) const;
    std::vector<double> Solve2(const std::vector<double> &b) const;

	friend std::ostream &operator<<(std::ostream &os, const Matrix &m);

	void Show() const;
    void LU(Matrix &L, Matrix &U, std::vector<double> &b) const;

    size_t eigenValuesJacobi(std::vector<double> &_eigenValuesvector, std::vector< std::vector<double> > &eigenVector, size_t itLimit = 1000) const;
    size_t eigenValuesScalar(double eigenValue, std::vector<double> &eigenVector, size_t itLimit = 1000) const;
    void eigenValuesDanilevski(std::vector<double> &eigenValuesvector, std::vector< std::vector<double> > &eigenVector) const;

    Matrix Trans() const;

    double sigma() const;
    double omega() const;
};