#include "SONE.h"
#include <fstream>
#include <algorithm>
#include "Matrix.h"

const double eps = 1E-5;
const long MAXIT = 1000;

extern std::ofstream cout;
using std::endl;

Item SONE::Norm(const vector<Item> &v) const
{
	Item r = Item(0);
	for (size_t i = 0; i < v.size(); i++)
		r += v[i] * v[i];

	return sqrt(r);
}

void SONE::ShowV(const vector<double> &v) const
{
	cout << "{";
	for (size_t i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << "}" << endl;
}

vector<Item> SONE::F(const vector<Item> &arg) const
{
	vector<Item> r;
	r.push_back(0.373 - sin(arg[1] - 1.938));
	r.push_back((-0.6 + cos(arg[0] + 0.061))/1.298);	

	return r;
}

vector<Item> operator-(const vector<Item> &v1, const vector<Item> &v2)
{
	size_t l = std::min(v1.size(), v2.size());
	vector<Item> r;
	r.resize(l);
	for (size_t i = 0; i < l; i++)
		r[i] = v1[i] - v2[i];

	return r;
}

vector<Item> SONE::f(const vector<Item> &arg) const
{
	vector<Item> r;
	r.push_back(-0.6 + cos(arg[0] + 0.061) - 1.298*arg[1]);
	r.push_back(-0.373 + sin(arg[1] - 1.938) + arg[0]);

	return r;
}

vector<Item> SONE::F2(const vector<Item> &arg) const
{
	vector<vector<Item>> Fder = { { cos(arg[0] + arg[1]) + 0.6, cos(arg[0] + arg[1]) }, { 2 * arg[0], 2 * arg[1] } };
	Matrix Fd(Fder);
	Fd = Fd.Inverse2();

	vector<Item> r1 = Fd*f2(arg);
	vector<Item> r = arg - r1;

	return r;
}

vector<Item> SONE::f2(const vector<Item> &arg) const
{
	vector<Item> r;
	r.push_back(-0.373 + sin(arg[0] + arg[1]) + 0.6*arg[0]);
	r.push_back(arg[0] * arg[0] + arg[1] * arg[1] - 1);

	return r;
}


vector<Item> SONE::SolveSimple()
{
	vector<Item> R = R_, residual = f(R), R0 = { 0, 0 };
	Item residualNorm = Norm(residual);
	long it = 1;

	while (residualNorm > eps || Norm(R - R0) > eps)
	{
		cout << "Iteration " << it << endl << "(x,y)=";
		ShowV(R);
		cout << "Residual vector=";
		ShowV(residual);
		cout << "Residual vector's norm="
			 << residualNorm << endl << endl;
		R0 = R;
		R = F(R);
		residual = f(R);
		residualNorm = Norm(residual);
		++it;
		if (it > MAXIT)
			throw std::exception("Too many iterations");
	} 

	cout << "Result in " << it << " iterations:" << endl << "(x,y)=";
	ShowV(R);
	cout << "Residual vector=";
	ShowV(residual);
	cout << "Residual vector's norm="
		<< residualNorm << endl;

	return R;
}

vector<Item> SONE::SolveNewton()
{
	vector<Item> R = R_, residual = f2(R), R0 = { 0, 0 };
	Item residualNorm = Norm(residual);
	long it = 1;

	while (residualNorm > eps || Norm(R - R0) > eps)
	{
		cout << "Iteration " << it << endl << "(x,y)=";
		ShowV(R);
		cout << "Residual vector=";
		ShowV(residual);
		cout << "Residual vector's norm="
			<< residualNorm << endl << endl;
		R0 = R;
		R = F2(R);
		residual = f2(R);
		residualNorm = Norm(residual);
		++it;
		if (it > MAXIT)
			throw std::exception("Too many iterations");
	}

	cout << "Result in " << it << " iterations:" << endl << "(x,y)=";
	ShowV(R);
	cout << "Residual vector=";
	ShowV(residual);
	cout << "Residual vector's norm="
		<< residualNorm << endl;

	return R;
}