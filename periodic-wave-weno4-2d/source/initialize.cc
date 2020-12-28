#include "../include/Weno432.h"


// Initialize the solution 

void Weno4_2D::initialize() {

    std::cout << "Initializing the solution" << std::endl;

    unsigned int N_gp = 4;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    Point<2> q_point;

    double V0;
    
    h_min = 10.0; 

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        if (c < no_cells_per_block) {

	        fv_values.reinit (cell);
		    V0 = 0.0;

			for (unsigned int i=0; i<fv_values.n_quadrature_points; ++i)
		  	  V0 += fv_values.JxW (i);

			double h = std::sqrt(V0);

            if (h < h_min) {            
                h_min = h; 
            }

            U(c) = 0.0;
            V(c) = 0.0; 

            for (unsigned int i = 0; i < N_gp*N_gp; i++) {

                q_point = fv_values.quadrature_point(i);

				U(c)   += (1./V0)*fv_values.JxW (i)*initial_condition_u(q_point);
                V(c)   += (1./V0)*fv_values.JxW (i)*initial_condition_v(q_point);  
            }
        }
    }
    
    std::cout << "Done!" << std::endl;  
    std::cout << "===========================" << std::endl;
} 
