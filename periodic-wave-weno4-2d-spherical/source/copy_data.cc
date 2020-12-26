#include "../include/Weno432.h"

// Copy periodic data 

void Weno4_2D::copy_data() {

    // Get iterators for cell	
	Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
	Triangulation<2>::active_cell_iterator endc = triangulation.end();
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        if (c >= no_cells_per_block && c < 2*no_cells_per_block) { // Bottom left 
            U(c)   = U(c - no_cells_per_block); 
        }
        
        if (c >= 2*no_cells_per_block && c < 3*no_cells_per_block) { // Bottom center 
            U(c)   = U(c - 2*no_cells_per_block); 
        }
        
        if (c >= 3*no_cells_per_block && c < 4*no_cells_per_block) {
            U(c)   = U(c - 3*no_cells_per_block); 
            
        }
        
        if (c >= 4*no_cells_per_block && c < 5*no_cells_per_block) {
            U(c)   = U(c - 4*no_cells_per_block); 
        }
        
        if (c >= 5*no_cells_per_block && c < 6*no_cells_per_block) {
            U(c)   = U(c - 5*no_cells_per_block); 
        }
        
        if (c >= 6*no_cells_per_block && c < 7*no_cells_per_block) {
            U(c)   = U(c - 6*no_cells_per_block); 
        }
        
        if (c >= 7*no_cells_per_block && c < 8*no_cells_per_block) {
            U(c)   = U(c - 7*no_cells_per_block); 
        }
        
        if (c >= 8*no_cells_per_block && c < 9*no_cells_per_block) {
            U(c)   = U(c - 8*no_cells_per_block); 
        }
    }
} 
