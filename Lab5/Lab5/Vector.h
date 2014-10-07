#pragma once

#include "Polynomial.h"
#include <vector>

vector<Item> operator+(const vector<Item> &v1, const vector<Item> &v2);
vector<Item> operator*(const vector<Item> v, double n);
vector<Item> operator*(double n, const vector<Item> v);
std::ostream &operator<<(std::ostream &os, const vector<Item> &v);