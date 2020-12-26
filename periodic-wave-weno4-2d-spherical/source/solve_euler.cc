#include "../include/Weno432.h"

void Weno4_2D::solve_euler()
{
	double time = 0.0;
	unsigned int count = 0;

	Vector<double> U_old(dof_handler.n_dofs());
	Vector<double> V_old(dof_handler.n_dofs());

	L_norm(time);

	while (time < finalTime)
	{
		compute_time_step_based_on_cfl_number(time);

		time += dt;
		
		std::cout << "time = " << time << ", Final time = " << finalTime<< std::endl;
		count++;

		compute_rhs();

		for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) 
		{
            if (c < no_cells_per_block) {
                U(c) = U(c) + dt*rhs1(c);
                V(c) = V(c) + dt*rhs2(c);
            }
		}

		 
		
		if (count%10 == 0){
			L_norm(time);
			output_results(count);
		} 
	}
}