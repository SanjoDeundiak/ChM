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
		 
	for (size_t i = 0; i < n; i++)
	{
		Item t;
		file >> t;
		b.push_back(t);
	}

	A->Solve2(b);

}