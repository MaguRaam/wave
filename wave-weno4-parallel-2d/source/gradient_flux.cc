#include "../include/gradient_flux.h"

double gradient_flux(double gLx, double gRx, double gLy, double gRy, double nx, double ny) {
	
	return 0.5*( (gLx + gRx)*nx  +  (gLy + gRy)*ny );	 
}