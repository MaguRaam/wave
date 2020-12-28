#include "../include/Weno432.h"

// Set the initial condtions u:
double Weno4_2D::initial_condition_u(Point<2> P)
{
    double x = P(0);
    double y = P(1);
    return sin(kx * x + ky * y);
}

double Weno4_2D::initial_condition_v(Point<2> P)
{
    double x = P(0);
    double y = P(1);
    return -omega * cos(kx * x + ky * y);
}

/** standing wave initial condition:
 * 
 * double x = P(0); double y = P(1); 
	int k = 1;
    return sin(M_PI * k * x) * sin(M_PI * k * y); 
 */
