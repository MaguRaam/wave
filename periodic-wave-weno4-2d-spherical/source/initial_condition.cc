#include "../include/Weno432.h"

// Set the initial condtions
double Weno4_2D::initial_condition(Point<2> P)
{
    double x = P(0);
    double y = P(1);
    double r = sqrt(x*x + y*y);
    return r < lambda ? (1.0/r)*cos(k*r) : 0.0;
}

/** standing wave initial condition:
 * 
 * double x = P(0); double y = P(1); 
	int k = 1;
    return sin(M_PI * k * x) * sin(M_PI * k * y); 
 */
