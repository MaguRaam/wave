#include "../include/Weno432.h"


// Perform the actual reconstruction 

void Weno4_2D::reconstruct() {
	
	// Variables for reconstruction of U

	unsigned int no_stencils = 10;
    unsigned int p = 4; 
    double epsilon = 1.0e-12; 
    double h; // Measure of cell size 
    
    double fourth_order_wt = 0.8;
	double third_order_wt = 0.03; 
	double second_order_wt = 0.0125; 
	double sum_gamma; 

	Vector<double> gamma(no_stencils);

	gamma(0) = fourth_order_wt; 
	gamma(1) = third_order_wt; 
	gamma(2) = third_order_wt;
	gamma(3) = third_order_wt;
	gamma(4) = third_order_wt;
	gamma(5) = third_order_wt;
	gamma(6) = second_order_wt;
	gamma(7) = second_order_wt;
	gamma(8) = second_order_wt;
	gamma(9) = second_order_wt;

    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;
    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;
    double nx, ny;   // Face normal vectors
     
    double u0;
    
	/* Fourth order stencil */ 
	Vector<double> d_u_4(4); 
    Vector<double> b_u_4; 
    Vector<double> u_coeff_4(9); 
    
	/* Third order centered stencil */ 
	Vector<double> d_u_3(4); 
	Vector<double> b_u_3; 
	Vector<double> u_coeff_3(5); 
    
	/* First third order one-sided stencil */ 
	Vector<double> b_u_31;         
	Vector<double> d_u_31;
	Vector<double> u_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_u_32;       
	Vector<double> d_u_32;
	Vector<double> u_coeff_32(5);
    
    /* Third third order one-sided stencil */
	Vector<double> b_u_33;          
	Vector<double> d_u_33;
	Vector<double> u_coeff_33(5);
    
    /* Fourth third order one-sided stencil */
    Vector<double> b_u_34; 
    Vector<double> d_u_34;
    Vector<double> u_coeff_34(5);   
	
 	Vector<double> u_coeff_21(2); 
	Vector<double> u_coeff_22(2); 
	Vector<double> u_coeff_23(2); 
	Vector<double> u_coeff_24(2);
	Vector<double> b_u2(2);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_U(no_stencils); Vector<double> w_U(no_stencils); double sum_U;
    

    // Iterate over all the cells 
    DoFHandler<2>::active_cell_iterator cell, neighbor;
    
    unsigned int index, ROWS, g_i; 

	unsigned int WW_index, NN_index, EE_index, SS_index;

	unsigned int neighbor_p1, neighbor_p2, neighbor_p3, neighbor_p4;
	
	for (unsigned int c = 0; c < n_relevant_cells; ++c) 
	{

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];

	    bool WW = false, EE = false, NN = false, SS = false;

        u0   =   U(g_i);
        h = std::sqrt(Cell[c].measure()); 

		if(!cell->face(0)->at_boundary()) {
	        neighbor = cell->neighbor(0);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p1 = local_neighbor_dof_indices[0];
		}

		if(!cell->face(3)->at_boundary()) {
	        neighbor = cell->neighbor(3);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p2 = local_neighbor_dof_indices[0];
		}

		if(!cell->face(1)->at_boundary()) {
	        neighbor = cell->neighbor(1);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p3 = local_neighbor_dof_indices[0];
		}

		if(!cell->face(2)->at_boundary()) {
	        neighbor = cell->neighbor(2);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p4 = local_neighbor_dof_indices[0];
		}

		if (cell_neighbor_neighbor_index[c][0].size() > 0) {
			WW = true;
			WW_index = cell_neighbor_neighbor_index[c][0][0] ;
		}

		if (cell_neighbor_neighbor_index[c][1].size() > 0) {
			EE = true;
			EE_index = cell_neighbor_neighbor_index[c][1][0] ;
		}

		if (cell_neighbor_neighbor_index[c][2].size() > 0) {
			SS = true;
			SS_index = cell_neighbor_neighbor_index[c][2][0] ;
		}

		if (cell_neighbor_neighbor_index[c][3].size() > 0) {
			NN = true;
			NN_index = cell_neighbor_neighbor_index[c][3][0] ;
		}
        
        if ( !(cell->at_boundary()) ) {
//			pcout<<"interior"<<std::endl;
			
            d_u_31.reinit(4);   d_u_32.reinit(4);   d_u_33.reinit(4);   
            d_u_34.reinit(4);
             
            
            // =====================================================================
            // r = 4 stencil 
            // =====================================================================
            
			d_u_4(0)  = (U(neighbor_p1) - u0);   // W neighbor 
			d_u_4(1)  = (U(neighbor_p2) - u0);   // N neighbor
			d_u_4(2)  = (U(neighbor_p3) - u0);   // E neighbor
			d_u_4(3)  = (U(neighbor_p4) - u0);   // S neighbor
			
            index = 0; 
            
            // Least Squares Part  
		
			ROWS = cell_diagonal_neighbor_index[c].size();
        			
			// neighbor of neighbors 

			if (WW) {
				ROWS++; 
			}

			if (EE) {
				ROWS++; 
			}

			if (SS) {
				ROWS++; 
			}

			if (NN) {
				ROWS++; 
			}
            
            b_u_4.reinit(ROWS);
			
			// vertex neighbors
			
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];

				b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
				index++; 
			}
			
			// neighbors of neighbors 
			
			if (WW) {

				local_neighbor_dof_indices[0] = WW_index;

				b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
				index++; 
			}
			
			if (NN) {

				local_neighbor_dof_indices[0] = NN_index;

				b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
				index++; 
			}
			
			if (EE) {

				local_neighbor_dof_indices[0] = EE_index;

				b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
				index++; 
			}
			
			if (SS) {

				local_neighbor_dof_indices[0] = SS_index;

				b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
				index++; 
			}

//			pcout<<"cls_R4 interior"<<std::endl;
			
			CLS_R4[c].solve(b_u_4, d_u_4, u_coeff_4);  
			
			// =====================================================================
            // r = 3 stencil (Centered stencil)
            // =====================================================================
			
			// constraint part (consists of face neighbours)
			
			d_u_3(0)  = (U(neighbor_p1) - u0);                // W neighbor 
			d_u_3(1)  = (U(neighbor_p2) - u0);                // N neighbor
			d_u_3(2)  = (U(neighbor_p3) - u0);                // E neighbor
			d_u_3(3)  = (U(neighbor_p4) - u0);                // S neighbor
			
					
			ROWS = cell_diagonal_neighbor_index[c].size();
			index = 0; 
            
            b_u_3.reinit(ROWS);
            
            /*TODO*/


			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];

				b_u_3(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
				 
				index++; 
			}

//			std::cout<<"cls_R3 interior"<<std::endl;			 
          
			CLS_R3[c].solve(b_u_3, d_u_3, u_coeff_3);

//			std::cout<<"cls_R3 interior end"<<std::endl;			 
			
			// =====================================================================
			// r = 3 stencil 1
			// =====================================================================
            
			if (is_admissible_R31[c]) {
            
				d_u_31(0)  = U(neighbor_p1) - u0;                // W neighbor 
				d_u_31(1)  = U(neighbor_p2) - u0;                // N neighbor
				d_u_31(2)  = U(neighbor_p4) - u0;                // S neighbor
				d_u_31(3)  = U(WW_index) - u0;               	 // WW neighbor
                
                ROWS = cell_neighbor_index[c][0].size(); 
				
				b_u_31.reinit(ROWS); 
                index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][0][d];

					b_u_31(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
					
					index++; 
				}

                CLS_R31[c].solve(b_u_31, d_u_31, u_coeff_31);
  
			}
            
            // If stencil is not available fallback to third order centered stencil 
            
            else {
				u_coeff_31 = u_coeff_3;
			}
            
            // =====================================================================
            // r = 3 stencil 2
            // =====================================================================
            
			if (is_admissible_R32[c]) {
            
				d_u_32(0)  = U(neighbor_p1) - u0;                // W neighbor 
				d_u_32(1)  = U(neighbor_p2) - u0;                // N neighbor
				d_u_32(2)  = U(neighbor_p3) - u0;                // E neighbor
				d_u_32(3)  = U(NN_index) -    u0;               // NN neighbor
                
                ROWS = cell_neighbor_index[c][3].size(); 

				b_u_32.reinit(ROWS); 
				index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][3][d];

					b_u_32(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
					index++; 
				}
//			std::cout<<"cls_R32 interior"<<std::endl;			 				
                CLS_R32[c].solve(b_u_32, d_u_32, u_coeff_32);
//			pcout<<"cls_R32 interior end"<<std::endl;			 				
			}
            // If stencil is not available fallback to third order centered stencil
            
            else {
				u_coeff_32 = u_coeff_3; 
            }  
            
            // =====================================================================
            // r = 3 stencil 3
            // =====================================================================
            
			if (is_admissible_R33[c]) {
            
				d_u_33(0)  = U(neighbor_p1) - u0;                // W neighbor 
				d_u_33(1)  = U(neighbor_p3) - u0;                // E neighbor
				d_u_33(2)  = U(neighbor_p4) - u0;                // S neighbor
				d_u_33(3)  = U(SS_index) -u0;               // SS neighbor
                
				

                ROWS = cell_neighbor_index[c][2].size(); 

				b_u_33.reinit(ROWS); 
                index = 0; 			

				for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][2][d];

					b_u_33(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
					index++; 
				}
//			std::cout<<"cls_R33 interior"<<std::endl;			 				
                CLS_R33[c].solve(b_u_33,   d_u_33,   u_coeff_33);
                 
//			pcout<<"cls_R33 interior end"<<std::endl;			 				
			}
            // If stencil is not available fallback to third order centered stencil
            
            else {
				u_coeff_33 = u_coeff_3; 
			} 
            
			// =====================================================================
            // r = 3 stencil 4
            // =====================================================================
            
			if (is_admissible_R34[c]) {
            
				d_u_34(0)  = U(neighbor_p2) - u0;                // N neighbor 
				d_u_34(1)  = U(neighbor_p3) - u0;                // E neighbor
				d_u_34(2)  = U(neighbor_p4) - u0;                // S neighbor
				d_u_34(3)  = U(EE_index) - u0;                // S neighbor
                

                ROWS = cell_neighbor_index[c][1].size(); 

                b_u_34.reinit(ROWS); 
                index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][1][d];

					b_u_34(index)   = U(local_neighbor_dof_indices[0] ) - u0;
					index++; 
				}
//			if (c == 241) pcout<<"cls_R34 interior rho"<<"\tc: "<<c<<std::endl<<b_rho_34<<std::endl<<d_rho_34<<std::endl;			 				
				CLS_R34[c].solve(b_u_34, d_u_34, u_coeff_34);			 				
            }
            
            // If stencil is not available fallback to third order centered stencil
            
            else {
				u_coeff_34 = u_coeff_3;
            } 
            
            
			// =====================================================================
			// r = 2 stencil 1
			// ===================================================================== 

			b_u2(0)  = (U(neighbor_p1) - u0);          // W neighbor 
			b_u2(1)  = (U(neighbor_p2) - u0);          // N neighbor
			LU_R21[c].solve(b_u2, u_coeff_21);
            
            // =====================================================================
            // r = 2 stencil 2 
            // ===================================================================== 

			b_u2(0) = (U(neighbor_p2) - u0);          // N neighbor
			b_u2(1) = (U(neighbor_p3) - u0);          // E neighbor
			LU_R22[c].solve(b_u2, u_coeff_22);

		
			// =====================================================================
			// r = 2 stencil 3
			// =====================================================================

			b_u2(0) = (U(neighbor_p3) - u0);          // E neighbor
			b_u2(1) = (U(neighbor_p4) - u0);          // S neighbor                                          
			LU_R23[c].solve(b_u2, u_coeff_23); 

			// =====================================================================
			// r = 2 stencil 4
			// =====================================================================

			b_u2(0) = (U(neighbor_p4) - u0);          // S neighbor
			b_u2(1) = (U(neighbor_p1) - u0);          // W neighbor
			LU_R24[c].solve(b_u2, u_coeff_24);
    
        } // End of interior cell loop 
        
        
        else {

            bool W_face = false, E_face = false, N_face = false, S_face = false;
            
            if (!(is_corner_cell[c])) {

//	   			pcout<<"boundary "<<std::endl;    
                
				if(cell->face(0)->at_boundary()) { W_face = true;}
				if(cell->face(1)->at_boundary()) { E_face = true;}
				if(cell->face(2)->at_boundary()) { S_face = true;}
				if(cell->face(3)->at_boundary()) { N_face = true;}
                
                // =====================================================================
                // r = 4 stencil (boundary)
                // =====================================================================
            
                Vector<double> d_u;  
                
                d_u.reinit(7);
                
                index = 0;
                
                if (!W_face) {
                    
                    d_u(index)   = U(neighbor_p1) - u0;
                    index++; 
                }
                
                if (!N_face) {
                    
                    d_u(index)   = U(neighbor_p2) - u0;
                    index++; 
                    
                }
                
                if (!E_face) {
                    
                    d_u(index)   = U(neighbor_p3) - u0;
                    index++; 
                    
                }
                
                if (!S_face) {
                    
                    d_u(index)   = U(neighbor_p4) - u0;
                    index++; 
                    
                }
                
                d_u(3)   = 0.0; d_u(4)   = 0.0; d_u(5)   = 0.0; d_u(6)   = 0.0;
                index = 0;
                
				// Least Squares Part  
			
				ROWS = cell_diagonal_neighbor_index[c].size(); 
				
				// neighbor of neighbors 
				
				if (WW) {
					ROWS++; 
				}
	
				if (EE) {
					ROWS++; 
				}
		
				if (SS) {
					ROWS++; 
				}

				if (NN) {
					ROWS++; 
				}
				
				b_u_4.reinit(ROWS);
				
				// vertex neighbors
				
				for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			

					local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];	

					b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
					index++; 
				}
								
				// neighbors of neighbors 
				
				if (WW) {

					local_neighbor_dof_indices[0] = WW_index;

					b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
					index++; 
				}
				
				if (NN) {

					local_neighbor_dof_indices[0] = NN_index;

					b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
					index++; 
				}
				
				if (EE) {

					local_neighbor_dof_indices[0] = EE_index;

					b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0;
					index++; 
				}
				
				if (SS) {

					local_neighbor_dof_indices[0] = SS_index;

					b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0;
					index++; 
				}
                
                
                CLS_R4[c].solve(b_u_4, d_u, u_coeff_4);

				// =====================================================================
                // r = 3 center stencil (boundary)
                // =====================================================================
                
                b_u_3.reinit(5);
                index = 0; 
                
                if (!W_face) {
                    
                    // P1 neighbor
                    b_u_3(index)    = (U(neighbor_p1) - u0);    
                    index++; 
                }
                
                if (!N_face) {

                    // P2 neighbor
                    
                    b_u_3(index)    = (U(neighbor_p2) - u0);
                    index++; 
                
                }
                
                if (!E_face) {

                    // P3 neighbor
                    
                    b_u_3(index)    = (U(neighbor_p3) - u0);
                    index++; 
                }
                
                if (!S_face) {

                    // P4 neighbor
                    
                    b_u_3(index)    = (U(neighbor_p4) - u0);
                    index++; 
                }
                
                
                // Transmissive boundary conditions 
                
				b_u_3(3) = 0.0; 
				b_u_3(4) = 0.0; 
                
            
				CLS_R3[c].solve(b_u_3, u_coeff_3);

				// =====================================================================
                // r = 3 stencil 1 (boundary)
                // =====================================================================
                
                if (is_admissible_R31[c]) {
                    
                    index = 0; 
                    if ( N_face || S_face ) {                    
                        d_u_31.reinit(5);   
                    }
            
                    else {
                        d_u_31.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_u_31(index)  = U(neighbor_p1) - u0;             
                        index++;
                    }
                    
                    if (!N_face) {
                    
                        d_u_31(index)  = U(neighbor_p2) - u0;
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_u_31(index)  = U(neighbor_p4) - u0; 
                        index++;
                    }
                    
                    if (WW) {

						local_neighbor_dof_indices[0] = WW_index;
                    
                        d_u_31(index)    = U(local_neighbor_dof_indices[0]) - u0; 
                        index++;
                    }
                    if ( S_face || N_face ) {                    
                        
                        // Transmissive boundary conditions 
                        
                        d_u_31(index)    = 0.0;   d_u_31(index+1)    = 0.0;        
                    }

	                ROWS = cell_neighbor_index[c][0].size(); 
				
					b_u_31.reinit(ROWS);
					index = 0; 
					
					// vertex neighbor of cell at face 0 
					for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][0][d];

						b_u_31(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
						index++; 
					}

					CLS_R31[c].solve(b_u_31, d_u_31, u_coeff_31);
					
				}
                
                else {
					u_coeff_31 = u_coeff_3; 
                }
                
				// =====================================================================
                // r = 3 stencil 2 (boundary)
                // =====================================================================
                
                if (is_admissible_R32[c]) {
                    
                    index = 0; 
                    if ( W_face || E_face ) {                    
                        d_u_32.reinit(5);  
                    }
            
                    else {
                        d_u_32.reinit(4); 
                    }
                    
                    if (!W_face) {
                    
                        d_u_32(index)  = U(neighbor_p1) - u0;
                        index++;
                    }
                    
                    if (!N_face) {
                    
                        d_u_32(index)  = U(neighbor_p2) - u0;                // N neighbor 
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_u_32(index)  = U(neighbor_p3) - u0;                // P3 neighbor 
                        index++;
                    }
                    
                    if (NN) {

						local_neighbor_dof_indices[0] = NN_index;
                     
                        d_u_32(index)    = U(local_neighbor_dof_indices[0] ) - u0; 
                        index++;
                    }
                    
                    if ( W_face || E_face ) {
                        
                        // Transmissive boundary conditions 
                        
                        d_u_32(index)    = 0.0;   d_u_32(index+1)    = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][3].size(); 

					b_u_32.reinit(ROWS);  

					index = 0; 		

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][3][d];

						b_u_32(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
						 
						index++; 
					}
                    
					CLS_R32[c].solve(b_u_32, d_u_32, u_coeff_32);
					
				}
                
                else {
					u_coeff_32 = u_coeff_3; 
					 
				}
                

                // =====================================================================
                // r = 3 stencil 3 (boundary)
                // =====================================================================
                
                if (is_admissible_R33[c]) {
                    
                    index = 0; 
                    if ( W_face || E_face ) {                    
                        d_u_33.reinit(5); 
                    }
            
                    else {
                        d_u_33.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_u_33(index)  = U(neighbor_p1) - u0;                // W neighbor 
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_u_33(index)  = U(neighbor_p3) - u0;                // E neighbor 
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_u_33(index)  = U(neighbor_p4) - u0;                // S neighbor 
                       
                        index++;
                    }
                    
                    if (SS) {

						local_neighbor_dof_indices[0] = SS_index;
                    
                        d_u_33(index)    = U(local_neighbor_dof_indices[0] ) - u0;              // SS neighbor 
                        
                        index++;
                    }
                    
                    if ( W_face || E_face ) {                    
                        
                        // Transmissive boundary conditions 
                        
                        d_u_33(index)    = 0.0;   d_u_33(index+1)    = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][2].size(); 

					b_u_33.reinit(ROWS);
					index = 0; 							

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][2][d];

						b_u_33(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
						
						index++; 
					}

					CLS_R33[c].solve(b_u_33, d_u_33, u_coeff_33);
					 
				}
                
                else {
					u_coeff_33 = u_coeff_3;
                }
                
                // =====================================================================
                // r = 3 stencil 4 (boundary)
                // =====================================================================
                
                if (is_admissible_R34[c]) {
                    
                    index = 0; 
                    if ( N_face || S_face ) {                                        
                        d_u_34.reinit(5);  
                    }
            
                    else {
                        d_u_34.reinit(4); 
                    }
                    
                    if (!N_face) {
                    
                        d_u_34(index)    = U(neighbor_p2) - u0;
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_u_34(index)    = U(neighbor_p3) - u0;              // E neighbor 
                        
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_u_34(index)    = U(neighbor_p4) - u0;              // S neighbor 
                        index++;
                    }
                    
                    if (EE) {

						local_neighbor_dof_indices[0] = EE_index;
                    
                        d_u_34(index)    = U(local_neighbor_dof_indices[0] ) - u0;              // S6 neighbor 
                        
                        index++;
                    }
                    
                    if ( N_face || S_face ) {                                        
                        
                        // Transmissive boundary conditions 
                        
                        d_u_34(index)    = 0.0;   d_u_34(index+1)    = 0.0; 
                    }
                    
	                ROWS = cell_neighbor_index[c][1].size(); 

					b_u_34.reinit(ROWS); 
					index = 0; 		

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][1][d];

						b_u_34(index)   = U(local_neighbor_dof_indices[0] ) - u0; 
						index++; 
					}                   
					
                    CLS_R34[c].solve(b_u_34, d_u_34, u_coeff_34);
                     
                }
                
                else {
					
					u_coeff_34 = u_coeff_3; 
                }
              
                // =====================================================================
                // r = 2 stencil 1 (boundary)
                // =====================================================================
                
                if (W_face) { // WEST boundary 
                
                    b_u2(0) = (U(neighbor_p2) - u0);   // P2 neighbor
                    b_u2(1) = 0.0; 
                    
                }
                
                else if (N_face) { // NORTH boundary 
                
                    b_u2(0) = (U(neighbor_p1) - u0);   // P2 neighbor
                    b_u2(1) = 0.0; 
                     
                }
                
                else {
                    
                    b_u2(0)  = (U(neighbor_p1) - u0);   // P1 neighbor 
                    b_u2(1)  = (U(neighbor_p2) - u0);   // P2 neighbor

                }
                
				LU_R21[c].solve(b_u2, u_coeff_21);

                // =====================================================================
                // r = 2 stencil 2 (boundary)
                // =====================================================================
                
                if (E_face) { // EAST boundary 
                
                    b_u2(0) = (U(neighbor_p2) - u0);   // P2 neighbor
                    b_u2(1) = 0.0; 
              
                }
                
                else if (N_face) { // NORTH boundary 
                
                    b_u2(0) = (U(neighbor_p3) - u0);   // P3 neighbor
                    b_u2(1) = 0.0; 
                }
                
                else {
                    
                    b_u2(0) = (U(neighbor_p2) - u0);   // P2 neighbor
                    b_u2(1) = (U(neighbor_p3) - u0);   // P3 neighbor
                    
                }
                
                LU_R22[c].solve(b_u2, u_coeff_22);
                // =====================================================================
                // r = 2 stencil 3 (boundary)
                // =====================================================================
                
                if (E_face) { // EAST boundary 
                
                    b_u2(0) = (U(neighbor_p4) - u0);   // P2 neighbor
                    b_u2(1) = 0.0; 
                    
                }
                
                else if (S_face) { // SOUTH boundary 
                
                    b_u2(0) = (U(neighbor_p3) - u0);   // P3 neighbor
                    b_u2(1) = 0.0; 
                    
                }
                
                else {
                
                    b_u2(0) = (U(neighbor_p3) - u0);   // P3 neighbor
                    b_u2(1) = (U(neighbor_p4) - u0);   // P4 
                }
                
                LU_R23[c].solve(b_u2, u_coeff_23);
                
                // =====================================================================
                // r = 2 stencil 4 (boundary)
                // =====================================================================
                
                if (W_face) { // WEST boundary 
                
                    b_u2(0) = (U(neighbor_p4) - u0);   // P4 neighbor
                    b_u2(1) = 0.0; 
                    
                }
                
                else if (S_face) { // SOUTH boundary 
                
                    b_u2(0) = (U(neighbor_p1) - u0);   // P1 neighbor
                    b_u2(1) = 0.0; 
                }
                
                else {
                    b_u2(0) = (U(neighbor_p4) - u0);   // P4 neighbor
                    b_u2(1) = (U(neighbor_p1) - u0);   // P1 neighbor
                }
                
                LU_R24[c].solve(b_u2, u_coeff_24);

				                    
            } // Non-corner boundary cell loop 
            
            else {
                
                // Corner boundary cells - reduce to first order 
                
                u_coeff_4 = 0.0;
				u_coeff_3 = 0.0;
				u_coeff_31 = 0.0;
				u_coeff_32 = 0.0;
                u_coeff_33 = 0.0;  
                u_coeff_34 = 0.0;  
				u_coeff_21 = 0.0;
				u_coeff_22 = 0.0;
				u_coeff_23 = 0.0;
				u_coeff_24 = 0.0;
				
				 
                
            } // Corner boundary cell loop 
            
        } // End of boundary cell loop 
        
        // Find the smoothness indicators 
        
		// =====================================================================
		// r = 4 Stencil 
		// =====================================================================
        
		IS_U(0)   = compute_fourth_order_smoothness_indicator(u_coeff_4, IS_constants[c], h);
		 
		// =====================================================================
		// r = 3
		// =====================================================================
        
		IS_U(1)   = compute_third_order_smoothness_indicator(u_coeff_3, IS_constants[c], h); 
		 
        
		// =====================================================================
		// r = 3  stencil 1
		// =====================================================================
        
		IS_U(2)   = compute_third_order_smoothness_indicator(u_coeff_31, IS_constants[c], h); 
		// =====================================================================
		// r = 3  stencil 2
		// =====================================================================
        
		IS_U(3)   = compute_third_order_smoothness_indicator(u_coeff_32, IS_constants[c], h); 
		
		// =====================================================================
		// r = 3  stencil 3
		// =====================================================================
        
		IS_U(4)   = compute_third_order_smoothness_indicator(u_coeff_33, IS_constants[c], h); 
        
		// =====================================================================
		// r = 3  stencil 4
		// =====================================================================
        
		IS_U(5)   = compute_third_order_smoothness_indicator(u_coeff_34, IS_constants[c], h); 
        
		// =====================================================================
		// r = 2  stencil 1
		// =====================================================================
		
		IS_U(6)   = compute_second_order_smoothness_indicator(u_coeff_21, IS_constants[c], h);
		// =====================================================================
		// r = 2  stencil 2
		// =====================================================================
		
		IS_U(7)   = compute_second_order_smoothness_indicator(u_coeff_22, IS_constants[c], h); 
		
		// =====================================================================
		// r = 2  stencil 3
		// =====================================================================
		
		IS_U(8)   = compute_second_order_smoothness_indicator(u_coeff_23, IS_constants[c], h); 
		// =====================================================================
		// r = 2  stencil 4
		// =====================================================================
		
		IS_U(9)   = compute_second_order_smoothness_indicator(u_coeff_24, IS_constants[c], h); 
		// Combine the polynomials 
		
		sum_U = 0.0; 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			
			w_U(j)   = gamma(j)/(std::pow((IS_U(j) + epsilon), p));
			
			sum_U += w_U(j);
			
			sum_gamma += gamma(j); 
		}
		
		// Normalize the weights 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			w_U(j) = w_U(j)/sum_U; 
			
			gamma(j) = gamma(j)/sum_gamma; 
		}
		
		// U
		
		coeffs_U[c](0) = u0; 
		
		coeffs_U[c](1) = (w_U(0)/gamma(0))*(u_coeff_4(0) - 
												
												gamma(1)*u_coeff_3(0)  - 
												gamma(2)*u_coeff_31(0) - 
												gamma(3)*u_coeff_32(0) - 
												gamma(4)*u_coeff_33(0) - 
												gamma(5)*u_coeff_34(0) - 
												gamma(6)*u_coeff_21(0) -
												gamma(7)*u_coeff_22(0) -
												gamma(8)*u_coeff_23(0) -
												gamma(9)*u_coeff_24(0) ) + 
												
												w_U(1)*u_coeff_3(0) +
												w_U(2)*u_coeff_31(0) + 
												w_U(3)*u_coeff_32(0) +
												w_U(4)*u_coeff_33(0) +
												w_U(5)*u_coeff_34(0) +
												w_U(6)*u_coeff_21(0) +
												w_U(7)*u_coeff_22(0) +
												w_U(8)*u_coeff_23(0) +
												w_U(9)*u_coeff_24(0) ;// u_x
		
		coeffs_U[c](2) = (w_U(0)/gamma(0))*(u_coeff_4(1) - 
												
												gamma(1)*u_coeff_3(1) - 
												gamma(2)*u_coeff_31(1) - 
												gamma(3)*u_coeff_32(1) -
												gamma(4)*u_coeff_33(1) -
												gamma(5)*u_coeff_34(1) - 
												gamma(6)*u_coeff_21(1) -
												gamma(7)*u_coeff_22(1) -
												gamma(8)*u_coeff_23(1) -
												gamma(9)*u_coeff_24(1) ) + 
												
												w_U(1)*u_coeff_3(1) +
												w_U(2)*u_coeff_31(1) +
												w_U(3)*u_coeff_32(1) +
												w_U(4)*u_coeff_33(1) +
												w_U(5)*u_coeff_34(1) +
												w_U(6)*u_coeff_21(1) +
												w_U(7)*u_coeff_22(1) +
												w_U(8)*u_coeff_23(1) +
												w_U(9)*u_coeff_24(1) ; // u_x

		coeffs_U[c](3) = (w_U(0)/gamma(0))*(u_coeff_4(2) - 
												
												gamma(1)*u_coeff_3(2) - 
												gamma(2)*u_coeff_31(2) - 
												gamma(3)*u_coeff_32(2) - 
												gamma(4)*u_coeff_33(2) - 
												gamma(5)*u_coeff_34(2) ) + 
												
												w_U(1)*u_coeff_3(2) +
												w_U(2)*u_coeff_31(2) +
												w_U(3)*u_coeff_32(2) + 
												w_U(4)*u_coeff_33(2) + 
												w_U(5)*u_coeff_34(2) ; // u_xx
												
		coeffs_U[c](4) = (w_U(0)/gamma(0))*(u_coeff_4(3) - 
												
												gamma(1)*u_coeff_3(3) - 
												gamma(2)*u_coeff_31(3) - 
												gamma(3)*u_coeff_32(3) - 
												gamma(4)*u_coeff_33(3) - 
												gamma(5)*u_coeff_34(3) ) + 
												
												w_U(1)*u_coeff_3(3) +
												w_U(2)*u_coeff_31(3) + 
												w_U(3)*u_coeff_32(3) + 
												w_U(4)*u_coeff_33(3) + 
												w_U(5)*u_coeff_34(3) ; // u_yy
		
		coeffs_U[c](5) = (w_U(0)/gamma(0))*(u_coeff_4(4) - 
												
												gamma(1)*u_coeff_3(4) - 
												gamma(2)*u_coeff_31(4) - 
												gamma(3)*u_coeff_32(4) - 
												gamma(4)*u_coeff_33(4) - 
												gamma(5)*u_coeff_34(4) ) + 
												
												w_U(1)*u_coeff_3(4) +
												w_U(2)*u_coeff_31(4) + 
												w_U(3)*u_coeff_32(4) + 
												w_U(4)*u_coeff_33(4) + 
												w_U(5)*u_coeff_34(4) ; // u_xy
		
		coeffs_U[c](6) = (w_U(0)/gamma(0))*(u_coeff_4(5)); // u_xxx

		coeffs_U[c](7) = (w_U(0)/gamma(0))*(u_coeff_4(6)); // u_yyy

		coeffs_U[c](8) = (w_U(0)/gamma(0))*(u_coeff_4(7)); // u_xxy

		coeffs_U[c](9) = (w_U(0)/gamma(0))*(u_coeff_4(8)); // u_xyy
		
 
        
    } // End of cell loop 
    
} // End of function 
