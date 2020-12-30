#include "../include/Weno432.h"

int main()
{

	std::cout.flags(std::ios::dec | std::ios::scientific);
	std::cout.precision(6);

	unsigned int cell = 128;

	double finalTime = 5.0;
	double cfl = 0.5;
	
	double c = 1.0;
	Weno4_2D test_problem(finalTime, cfl, cell, c);
	test_problem.run();

	return 0;
}
