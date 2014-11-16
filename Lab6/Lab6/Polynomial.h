#pragma once

#include <vector>
#include <iostream>

using std::vector;

typedef double Item;

class Polynomial
{
	vector<Item> vec; // from x^0 to x^n
public:
	// Constructors
	explicit Polynomial(const vector<Item> &A);

	Item operator()(Item x) const; // Gives value in point x
	Polynomial Derivative() const;

	// Arithmetic operators
	Polynomial operator*(Item n) const;
	friend Polynomial operator*(Item n, const Polynomial p) { return p*n; }
	Polynomial operator*(const Polynomial &other) const;
	Polynomial Polynomial::operator+(const Polynomial &other) const;	
	vector<Item> C() const { return vec;  };
    Item Polynomial::GetRootBin(double left, double right, double acc) const;

	friend std::ostream &operator<<(std::ostream &os, const Polynomial &p);
};