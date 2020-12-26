#include "../include/Weno432.h"

void Weno4_2D::solve_ssprk33()
{
	double time = 0.0;
	unsigned int count = 0;

	Vector<double> U_old(dof_handler.n_dofs());
	Vector<double> V_old(dof_handler.n_dofs());

	L_norm(time);

	while (time < finalTime)
	{
		compute_time_step_based_on_cfl_number(time);

		if(count == 0){
		    copy_data();
		}


		time += dt;
		
		std::cout << "time = " << time << ", Final time = " << finalTime<< std::endl;
		count++;

		U_old = U;
		V_old = V;

		// SSPRK Step 1

		compute_rhs();

		for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) 
		{
            if (c < no_cells_per_block) {
                U(c) = U(c) + dt*rhs1(c);
                V(c) = V(c) + dt*rhs2(c);
            }
		}

		// SSPRK Step 2
		compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
            if (c < no_cells_per_block) {
                U(c)   = (3./4.)*U_old(c)   + (1./4.)*U(c)   + (1./4.)*dt*rhs1(c);
                V(c)   = (3./4.)*V_old(c)   + (1./4.)*V(c)   + (1./4.)*dt*rhs2(c);
            }
        }

        // SSPRK Step 3
		compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
            if (c < no_cells_per_block) {
                U(c)   = (1./3.)*U_old(c)   + (2./3.)*U(c)   + (2./3.)*dt*rhs1(c);
                V(c)   = (1./3.)*V_old(c)   + (2./3.)*V(c)   + (2./3.)*dt*rhs2(c);
            }
		}

		L_norm(time);
		if (count % 10 == 0)  output_results(count);
 
	}//end of time loop

	 
}