#include "../include/Weno432.h"


double exact_solution(Point<2> P, double t)
{
	int k = 1;
	double omega = sqrt(2*M_PI*M_PI*k*k);

	double x = P(0); double y = P(1);
	return cos(omega * t) * sin(M_PI * k * x) * sin(M_PI * k * y);
}