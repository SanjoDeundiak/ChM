#pragma once

#include <vector>
#include <string>
#include <fstream>

using std::vector;
using std::string;
using std::ifstream;

typedef double Item;

class Matrix
{
	vector <vector <Item>> A;
	size_t n;

	vector<Item> SolveL(const vector<Item> &b) const;
	vector<Item> SolveU(const vector<Item> &b) const;
public:
	Matrix(size_t dimension, Item fill);
	Matrix(size_t dimension);
	Matrix(size_t dimension, ifstream &file);	

	Item Determinant() const;

	Matrix Inverse();

	Matrix operator*(const Matrix &other) const;
	vector<Item> operator*(const vector<Item> &vec) const;

	vector<Item> Solve(const vector<Item> &b) const;

	void Show() const;
	void LU(Matrix &L, Matrix &U) const;
};


