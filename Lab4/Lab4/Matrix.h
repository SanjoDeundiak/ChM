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

public:
	Matrix(const vector<vector<Item>> &v);	
	Matrix(size_t n, Item fill);

	Matrix Inverse2() const;

	Matrix operator*(const Matrix &other) const;
	vector<Item> operator*(const vector<Item> &vec) const;		
};


