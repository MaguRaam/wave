#include "../include/Weno432.h"


//  Evaluate time step using the CFL condition 

void Weno4_2D::compute_time_step_based_on_cfl_number(double time) {

	dt = 0.5*(cfl*h_min)/C;

	if((time + dt)>finalTime) {
		dt = finalTime- time;
	}
} 