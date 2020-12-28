#include "../include/Weno432.h"

int main()
{

	std::cout.flags(std::ios::dec | std::ios::scientific);
	std::cout.precision(6);

	unsigned int cell = 128;

	double finalTime = 5.0;
	double cfl = 0.5;

	//plane wave parameters:
	const double lambda(1.0);							 //wave length
	const double k(2.0 * M_PI / lambda);				 //wave number
	const double alpha(0.0);							 //wave propagation angle in radians
	const double kx(k * cos(alpha)), ky(k * sin(alpha)); //wave number x and y components:
	const double T(1.0), omega(2.0 * M_PI / T);			 //time period and angular frequency:

	//dispersion relation to compute wave speed:
	const double c(omega / k);
	std::cout<<"wave speed = "<<c<<std::endl;

	Weno4_2D test_problem(finalTime, cfl, cell, kx, ky, omega, c);
	test_problem.run();

	return 0;
}
