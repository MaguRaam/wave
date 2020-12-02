#include "../include/Weno432.h"

void Weno4_2D::L_norm(const double time) {

	DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    QGauss<2> quadrature_formula(4);

    const UpdateFlags update_flags =   update_values | update_quadrature_points | update_JxW_values;
	FEValues<2> fv_values (fv, quadrature_formula, update_flags);

	const unsigned int n_q_points = fv_values.n_quadrature_points;
	 

	for (unsigned int c = 0; cell != endc; ++cell, ++c) {
		
		if (c < no_cells_per_block) {

			fv_values.reinit(cell);

			double avg_exact_value = 0.0;
			double volume = 0.0;

			for (unsigned int i=0; i<n_q_points; ++i)
				volume += fv_values.JxW (i);

			for (unsigned int i=0; i<n_q_points; ++i) {
				
				Point<2> quad_point;
				quad_point = fv_values.quadrature_point(i);
	
				avg_exact_value += (1.0/volume)*exact_solution(quad_point,time)*fv_values.JxW (i);
			}

			local_difference(c) = avg_exact_value - U[c];

		}
	}

	std::ofstream fout_convergence ; 
	fout_convergence.flags( std::ios::dec | std::ios::scientific ) ; 
	fout_convergence.precision(6) ;

    const std::string filename = "../error/norm.dat";
	fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

   	fout_convergence << time 
   					<<"\t"<<local_difference.linfty_norm()	
   					<<"\t"<<local_difference.l2_norm()/std::sqrt(no_cells_per_block)
   					<< std::endl;
	fout_convergence.close();

}