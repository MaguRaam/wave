#include "../include/Weno432.h"


// Put everything together - solve the actual problem  

void Weno4_2D::run() {
    
    auto start = std::chrono::system_clock::now();
        
    make_grid(); 
	setup_system();
	compute_cell_properties();
    compute_weno_polynomial_constants(); 
	initialize();
	precompute_matrices(); 
    compute_IS_constants(); 
	solve_ssprk33(); 
   	 
    
    auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "Time taken = " << elapsed_seconds.count() << std::endl;
} 
