#pragma once

#include "Polynomial.h"
#include "Spline.h"
#include <vector>

using std::vector;

class Interpolation
{
	Item OriginalF(Item x) const;
	vector<Item> X_;
	vector<Item> Y_;
	Polynomial Newton1_(const vector<Item> &X, const vector<Item> &Y) const;
	Polynomial Lagrange_(const vector<Item> &X, const vector<Item> &Y) const;

public:
	Interpolation();
	Interpolation(Item a, Item b, int n);
	Polynomial Lagrange() const;
	Polynomial Newton1() const;
	Polynomial Newton2() const;
	Item Left() const { return X_[0]; };
	Item Right() const { return X_[X_.size() - 1]; };
	size_t n() const { return X_.size(); };
	Spline Spline2() const;
};