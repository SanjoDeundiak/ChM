#include "SOLE.h"
#include <iostream>
#include <string>
#include <fstream>

const double eps = 1E-5;

using std::cin;
//using std::cout;
using std::endl;
using std::ifstream;

std::ofstream cout("output.txt");

const string INFILE = "input.txt";

void SOLE::start()
{
	ifstream file(INFILE);
	size_t n;
	if (!file.is_open())
		throw std::runtime_error("Can't open file");

	file >> n;
	A = std::make_shared<Matrix>(n, file);
	cout << "Input Matrix" << endl;
	A->Show();

	Item det = A->Determinant();
	if (abs(det) < eps)
		throw std::exception("Singular matrix");

	cout << "Determinant=" << det << endl;

	cout << "A^(-1): " << endl;
	A->Inverse().Show();

	cout << "A*A^(-1):" << endl;
	(A->Inverse()*(*A)).Show();
	 
	for (size_t i = 0; i < n; i++)
	{
		Item t;
		file >> t;
		b.push_back(t);
	}

	vector<Item> r = A->Solve(b);
	cout << "Result={";
	for (size_t i = 0; i < n; i++)
		cout << r[i] << ' ';
	cout << '}' << endl;

	cout << "Ax-b={";
	r = *A * r;
	for (size_t i = 0; i < n; i++)
		cout << r[i] - b[i] << ' ';
	cout << '}' << endl;

}