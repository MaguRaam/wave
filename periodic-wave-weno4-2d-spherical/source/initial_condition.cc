#include "../include/Weno432.h"

// Set the initial condtions
double Weno4_2D::initial_condition(Point<2> P)
{
    double x = P(0);
    double y = P(1);
    
    PetscReal xs(0.0), ys(0.0); //source location:
    PetscReal w(0.125);      //half width:
   auto r = ((x - xs) * (x - xs) + (y - ys) * (y - ys));
    return std::exp(-std::log(2.0) * (r * r) / (w * w));
}

/** standing wave initial condition:
 * 
 * double x = P(0); double y = P(1); 
	int k = 1;
    return sin(M_PI * k * x) * sin(M_PI * k * y); 
 */
