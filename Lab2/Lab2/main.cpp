#include "SOLE.h"
#include <iostream>


using std::cin;
using std::cout;

int main()
{
	SOLE sole;
	try
	{
		sole.start();
	}
	catch (std::exception e)
	{
		cout << e.what();
	}

	cout << "DONE";
	cin.get();
	cin.get();
}