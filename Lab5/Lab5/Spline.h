#pragma once
#include"Polynomial.h"
#include <vector>
using std::vector;

class Spline
{
	vector<Polynomial> P;
	vector<Item> X;
	
public:
	Spline(Item x);
	void Add(const Polynomial &p, Item x);

	Item operator()(Item x) const; // Gives value in point x
	friend std::ostream &operator<<(std::ostream &os, const Spline &p);
};