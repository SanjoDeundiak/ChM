#pragma once
#include "Matrix.h"
#include <memory>

using std::shared_ptr;

class SOLE
{
	shared_ptr<Matrix> A;
	vector<Item> b;
public:
	void start();
};