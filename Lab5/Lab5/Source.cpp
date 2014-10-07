#include <fstream>
#include <array>
#include "Matrix.h"

using namespace std;

int main()
{
	ofstream out("output_shad.txt");
	vector<int> c(9,0);		
	vector<vector<double>> c1;
	vector<double> temp;

	for (int i = 0; i < 1953125; i++)
	{
		c1.clear();
		temp.clear();
		for (int j = 0; j < 25; j++)
		{
			temp.push_back(c[j]);
			if ((j + 1) % 5 == 0)
			{
				c1.push_back(temp);
				temp.clear();
			}
		}
		Matrix a(c1);
		out << a << endl << endl;
		/*double d;
		try
		{
			d = a.Determinant();
		}
		catch (...)
		{
			d = 0;
		}
		out << d << endl << endl;*/
		int j = 24;
		for (; c[j] == 1 && j >= 0; c[j] = 0, j--);
		if (j >= 0)
			c[j] = 1;
	}
	
}