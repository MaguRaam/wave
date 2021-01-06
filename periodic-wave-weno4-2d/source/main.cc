#include "../include/Weno432.h"

int main()
{

	std::cout.flags(std::ios::dec | std::ios::scientific);
	std::cout.precision(6);

	unsigned int cell = 256;

	double cfl = 0.5;

	//plane wave parameters:
	const double kx(2.0 * M_PI / 4.0), ky(2.0 * M_PI / 4.0); //wave number x and y components:
	const double T(1.0), omega(2.0 * M_PI / T);				 //time period and angular frequency:

	//dispersion relation to compute wave speed:
	const double c = omega / sqrt(kx * kx + ky * ky);
	std::cout << "wave speed = " << c << std::endl;

	double finalTime = T;
	Weno4_2D test_problem(finalTime, cfl, cell, kx, ky, omega, c);
	test_problem.run();

	return 0;
}
