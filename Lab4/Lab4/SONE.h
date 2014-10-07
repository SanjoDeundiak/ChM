#pragma once

#include <vector>

using std::vector;

typedef double Item;

class SONE
{
	const size_t n;
	vector<Item> F(const vector<Item> &arg) const;
	vector<Item> f(const vector<Item> &arg) const;

	vector<Item> F2(const vector<Item> &arg) const;
	vector<Item> f2(const vector<Item> &arg) const;

	void ShowV(const vector<double> &v) const;

	Item Norm(const vector<Item> &v) const;

	vector<Item> R_;
	SONE() : n(0) {};
public:
	SONE(Item x0, Item y0) : n(2)
	{		
		R_.push_back(x0);
		R_.push_back(y0);		
	}
	vector<Item> SolveSimple();
	vector<Item> SolveNewton();

};