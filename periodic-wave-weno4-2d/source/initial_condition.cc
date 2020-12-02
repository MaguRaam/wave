#include "../include/Weno432.h"

// Set the initial condtions  
double initial_condition(Point <2> P) {
	
	double x = P(0); double y = P(1); 
	int k = 1;
    return sin(M_PI * k * x) * sin(M_PI * k * y);
}
