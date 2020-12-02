#include "../include/Weno432.h"

// Setup the system - allocate the necessary memory 

void Weno4_2D::allocate_memory() {

    pcout << "allocate memory" << std::endl;

	U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	
	local_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	
    coeffs_U.resize(n_relevant_cells);
	
	Cell.resize(n_store_cell);
    
    WENO_poly_consts.resize(n_relevant_cells); 
    IS_constants.resize(n_relevant_cells);
    is_corner_cell.resize(n_relevant_cells);
    
    CLS_R4.resize(n_relevant_cells);
    CLS_R3.resize(n_relevant_cells);
    
    is_admissible_R31.resize(n_relevant_cells);
    CLS_R31.resize(n_relevant_cells);
    
	is_admissible_R32.resize(n_relevant_cells);
    CLS_R32.resize(n_relevant_cells);
    
    is_admissible_R33.resize(n_relevant_cells);
    CLS_R33.resize(n_relevant_cells);
    
	is_admissible_R34.resize(n_relevant_cells);
    CLS_R34.resize(n_relevant_cells);
    
	LU_R21.resize(n_relevant_cells);
    LU_R22.resize(n_relevant_cells);
    LU_R23.resize(n_relevant_cells);
    LU_R24.resize(n_relevant_cells);

	rhs1.reinit(n_locally_cells);
	rhs2.reinit(n_locally_cells);
	 
    
    for (unsigned int i = 0; i < n_relevant_cells; i++) {
        coeffs_U[i].reinit(10);
        WENO_poly_consts[i].reinit(9); 
        IS_constants[i].reinit(15);
    }

	Utilities::System::MemoryStats stat;
	Utilities::System::get_memory_stats(stat);
	pcout<<stat.VmRSS/std::pow(2,20)<<std::endl;
	pcout<<"total (on 24 core) : "<<24.0*stat.VmRSS/std::pow(2,20)<<std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "log.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Memory consumption in Allocate memory per node in GB = " << stat.VmRSS/std::pow(2,20) << std::endl;
        fout_convergence.close();
	}

}
