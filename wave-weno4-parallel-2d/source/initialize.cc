#include "../include/Weno432.h"


// Initialize the solution 

void Weno4_2D::initialize() 
{

    unsigned int N_gp = 4;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator cell;

    Point<2> q_point;

    double V0;

	double h_old,h_min_local; 
	h_min_local = 1e6; 

	unsigned int g_i;

	for (unsigned int c = 0; c < n_locally_cells; ++c) 
	{

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c]; 

		fv_values.reinit(cell);

		V0 = cell->measure(); 
	
		h_old = std::sqrt(V0);

		if (h_old < h_min_local) {
			h_min_local = h_old; 
		}

		double u = 0.0, v = 0.0;

		for (unsigned int i = 0; i < N_gp*N_gp; i++) {

			q_point = fv_values.quadrature_point(i);
			u += (1./V0)*fv_values.JxW (i)*initial_condition(q_point);
		}

		local_U(g_i) = u;
		local_V(g_i) = v;

    }

	local_U.compress(VectorOperation::insert);
	local_V.compress(VectorOperation::insert);

	U = local_U;
	V = local_V;
 
	h_min = Utilities::MPI::min (h_min_local, MPI_COMM_WORLD);
    pcout << "h_min: " <<h_min<< std::endl;  

	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
	    std::ofstream fout_convergence ;
	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
	    fout_convergence.precision(7) ;

	    const std::string filename = "log.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

		fout_convergence << "h_min: " <<h_min<< std::endl;  
	    fout_convergence.close();
	}

} 
