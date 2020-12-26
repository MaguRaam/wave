#include "../include/Weno432.h"


double Weno4_2D::exact_solution(Point<2> P, double t)
{
	double x = P(0); double y = P(1);
	return sin(kx*x + ky*y - omega*t);
}