#include "../include/Weno432.h"

int main()
{

	std::cout.flags(std::ios::dec | std::ios::scientific);
	std::cout.precision(6);

	unsigned int cell = 64;

	double finalTime = 5.0;
	double cfl = 0.5;

	//plane wave parameters:
	const double lambda(0.3);							 		//wave length
	const double k(2.0 * M_PI / lambda);				   //wave number
	const double T(1.0), omega(2.0 * M_PI / T);			 //time period and angular frequency:

	//dispersion relation to compute wave speed:
	const double c(omega / k);
	std::cout<<"wave speed = "<<c<<std::endl;

	Weno4_2D test_problem(finalTime, cfl, cell, k, lambda, omega, c);
	test_problem.run();

	return 0;
}
