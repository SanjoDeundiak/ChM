#include "SONE.h"
#include <fstream>
#include <iostream>

std::ofstream cout("output.txt");

int main()
{
	SONE S(1.2, 0), S1(0, 1), S2(1,-0.7);
	try
	{
		S2.SolveNewton();
	}
	catch (std::exception &e)
	{
		std::cout << e.what() << std::endl;
	}

	std::cin.get();
}