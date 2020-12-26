#include "../include/Weno432.h"


// Setup the system - allocate the necessary memory 

void Weno4_2D::setup_system() {

    std::cout << "Setting up system" << std::endl;

	dof_handler.distribute_dofs (fv);

	U.reinit(dof_handler.n_dofs());

    V.reinit(dof_handler.n_dofs());

	rhs1.reinit(no_cells_per_block);

    rhs2.reinit(no_cells_per_block);

 

    Cell.resize(no_cells_per_block); 

	local_difference.reinit(dof_handler.n_dofs());
    Uexact.reinit(dof_handler.n_dofs());	

    coeffs_U.resize(no_cells_per_block);
    
    WENO_poly_consts.resize(no_cells_per_block); 
    IS_constants.resize(no_cells_per_block);
    
    CLS_R4.resize(no_cells_per_block);

    CLS_R3.resize(no_cells_per_block);
    
    is_admissible_R31.resize(no_cells_per_block);
    CLS_R31.resize(no_cells_per_block);
    
	is_admissible_R32.resize(no_cells_per_block);
    CLS_R32.resize(no_cells_per_block);
    
    is_admissible_R33.resize(no_cells_per_block);
    CLS_R33.resize(no_cells_per_block);
    
	is_admissible_R34.resize(no_cells_per_block);
    CLS_R34.resize(no_cells_per_block);
    
	LU_R21.resize(no_cells_per_block);
    LU_R22.resize(no_cells_per_block);
    LU_R23.resize(no_cells_per_block);
    LU_R24.resize(no_cells_per_block);
    
    for (unsigned int i = 0; i < no_cells_per_block; i++) {
        coeffs_U[i].reinit(10);
        WENO_poly_consts[i].reinit(9); 
        IS_constants[i].reinit(15);
    }

    std::cout << "Done!" << std::endl;  
    std::cout << "===========================" << std::endl;
}
