#include "../include/Weno432.h"


int main() {

	std::cout.flags( std::ios::dec | std::ios::scientific ) ; 
	std::cout.precision(6) ;

	unsigned int cell = 64;

	double finalTime = 1.0;
	double cfl = 0.45;
	Weno4_2D test_problem(finalTime, cfl, cell);
	test_problem.run();

    return 0;
} 
