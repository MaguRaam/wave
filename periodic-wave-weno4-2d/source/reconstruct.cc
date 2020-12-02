#include "../include/Weno432.h"


// Perform the actual reconstruction 

void Weno4_2D::reconstruct() {
    
    // Variables for reconstruction of U
    
    unsigned int no_stencils = 10;
    unsigned int p = 4; 
    double epsilon = 1.0e-12; 
    double h; // Measure of cell size 
    
    double fourth_order_wt = 8.00;
	double third_order_wt = 0.4; 
	double second_order_wt = 0.0; 
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
    Vector<double> IS_u(no_stencils); Vector<double> w_u(no_stencils); double sum_u;
    
	// Iterate over all the cells 
    
    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    
    unsigned int ROWS, index; 
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        if (c < no_cells_per_block) {
        
            u0   =   U(cell->active_cell_index());
	        h = std::sqrt(Cell[c].measure()); 

            d_u_31.reinit(4); d_u_32.reinit(4); d_u_33.reinit(4); d_u_34.reinit(4);
            
            // =====================================================================
            // r = 4 stencil 
            // =====================================================================
            
            d_u_4(0)  = (U(cell->neighbor_index(0)) - u0);   // P1 neighbor 
            d_u_4(1)  = (U(cell->neighbor_index(3)) - u0);   // P2 neighbor
            d_u_4(2)  = (U(cell->neighbor_index(1)) - u0);   // P3 neighbor
            d_u_4(3)  = (U(cell->neighbor_index(2)) - u0);   // P4 neighbor
		
            ROWS = 0; index = 0; 
            
            // Least Squares Part  
        
            
            if (cell->neighbor(0)->neighbor_index(0) != -1) {
                ROWS++; // S1 Cell 
            }
            
            if (cell->neighbor(0)->neighbor_index(3) != -1) {
                ROWS++; // S2 cell
            }
            
            if (cell->neighbor(3)->neighbor_index(3) != -1) {
                ROWS++; // S3 Cell 
            }
            
            if (cell->neighbor(3)->neighbor_index(1) != -1) {
                ROWS++; // S4 Cell 
            }
            
            if (cell->neighbor(1)->neighbor_index(1) != -1) {
                ROWS++; // S5 Cell 
            }
            
            if (cell->neighbor(1)->neighbor_index(2) != -1) {
                ROWS++; // S6 Cell 
            }
            
            if (cell->neighbor(2)->neighbor_index(2) != -1) {
                ROWS++; // S7 Cell 
            }
            
            if (cell->neighbor(2)->neighbor_index(0) != -1) {
                ROWS++; // S8 Cell 
            }

            
            b_u_4.reinit(ROWS);
            
            if (cell->neighbor(0)->neighbor_index(0) != -1) {
                b_u_4(index)   = U(cell->neighbor(0)->neighbor_index(0)) - u0; 
                index++; 
            }
            
            if (cell->neighbor(0)->neighbor_index(3) != -1) {
                b_u_4(index)   = U(cell->neighbor(0)->neighbor_index(3)) - u0; 
                index++; 
            }
            
            if (cell->neighbor(3)->neighbor_index(3) != -1) {
                b_u_4(index)   = U(cell->neighbor(3)->neighbor_index(3)) - u0; 
                index++; 
            }
            
            if (cell->neighbor(3)->neighbor_index(1) != -1) {
                b_u_4(index)   = U(cell->neighbor(3)->neighbor_index(1)) - u0; 
                index++; 
            }
            
            if (cell->neighbor(1)->neighbor_index(1) != -1) {
                b_u_4(index)   = U(cell->neighbor(1)->neighbor_index(1)) - u0; 
                index++; 
            }
            
            if (cell->neighbor(1)->neighbor_index(2) != -1) {
                b_u_4(index)   = U(cell->neighbor(1)->neighbor_index(2)) - u0; 
                index++; 
            }
            
            if (cell->neighbor(2)->neighbor_index(2) != -1) {
                b_u_4(index)   = U(cell->neighbor(2)->neighbor_index(2)) - u0; 
                index++; 
            }
            
            if (cell->neighbor(2)->neighbor_index(0) != -1) {
                b_u_4(index)   = U(cell->neighbor(2)->neighbor_index(0)) - u0; 
            }
            
            CLS_R4[c].solve(b_u_4, d_u_4, u_coeff_4); 

			// =====================================================================
            // r = 3 stencil (Centered stencil)
            // =====================================================================
            
            ROWS = 0; index = 0; 
            
            if (cell->neighbor(0)->neighbor_index(3) != -1) {
                ROWS++; 
            }

            if (cell->neighbor(3)->neighbor_index(1) != -1) {
                ROWS++; 
            }
            
            if (cell->neighbor(1)->neighbor_index(2) != -1) {
                ROWS++; 
            }

            if (cell->neighbor(2)->neighbor_index(0) != -1) {
                ROWS++; 
            }
            
            b_u_3.reinit(ROWS);

            d_u_3(0)  = (U(cell->neighbor_index(0)) - u0);                // P1 neighbor 
            d_u_3(1)  = (U(cell->neighbor_index(3)) - u0);                // P2 neighbor
            d_u_3(2)  = (U(cell->neighbor_index(1)) - u0);                // P3 neighbor
            d_u_3(3)  = (U(cell->neighbor_index(2)) - u0);                // P4 neighbor
            
            if (cell->neighbor(0)->neighbor_index(3) != -1) {
                
                // S2 cell 
                
                b_u_3(index) = (U(cell->neighbor(0)->neighbor_index(3)) - u0);   
                index++; 
            }

            if (cell->neighbor(3)->neighbor_index(1) != -1) {
                
                // S4 cell 
                
                b_u_3(index) = (U(cell->neighbor(3)->neighbor_index(1)) - u0);    
                index++; 
            
            }
            
            if (cell->neighbor(1)->neighbor_index(2) != -1) {
                
                // S6 cell 
                
                b_u_3(index) = (U(cell->neighbor(1)->neighbor_index(2)) - u0);     
                index++;  
            }

            if (cell->neighbor(2)->neighbor_index(0) != -1) {
                
                // S8 cell 
                
                b_u_3(index) = (U(cell->neighbor(2)->neighbor_index(0)) - u0);   
  
            }
            
			CLS_R3[c].solve(b_u_3, d_u_3, u_coeff_3);
			
			// =====================================================================
			// r = 3 stencil 1
			// =====================================================================
            
			if (is_admissible_R31[c]) {
            
				d_u_31(0)  = U(cell->neighbor_index(0)) - u0;                // P1 neighbor 
				d_u_31(1)  = U(cell->neighbor_index(3)) - u0;                // P2 neighbor
				d_u_31(2)  = U(cell->neighbor_index(2)) - u0;                // P4 neighbor
				d_u_31(3)  = U(cell->neighbor(0)->neighbor_index(0)) - u0;   // S1 neighbor
                
				if (cell->neighbor(0)->neighbor_index(3) != -1  &&
					cell->neighbor(2)->neighbor_index(0) != -1) {
					ROWS = 2; 
				}

                else {
                    ROWS = 1; 
                }

                b_u_31.reinit(ROWS);

                index = 0; 

                if (cell->neighbor(0)->neighbor_index(3) != -1) { // S2 cell
                    b_u_31(index)   =  U(cell->neighbor(0)->neighbor_index(3)) - u0;  
                    index++; 
                }

                if (cell->neighbor(2)->neighbor_index(0) != -1) { // S8 cell
                    b_u_31(index)   =  U(cell->neighbor(2)->neighbor_index(0)) - u0;
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
            
                // Form constraint RHS

				d_u_32(0)  = U(cell->neighbor_index(0)) - u0;                // P1 neighbor 
				d_u_32(1)  = U(cell->neighbor_index(3)) - u0;                // P2 neighbor
				d_u_32(2)  = U(cell->neighbor_index(1)) - u0;                // P3 neighbor
				d_u_32(3)  = U(cell->neighbor(3)->neighbor_index(3)) - u0;   // S3 neighbor

                if (cell->neighbor(0)->neighbor_index(3) != -1  &&
                    cell->neighbor(3)->neighbor_index(1) != -1) {
					ROWS = 2; 
                }

                else {
					ROWS = 1; 
                }

				b_u_32.reinit(ROWS);

                index = 0; 

                if (cell->neighbor(0)->neighbor_index(3) != -1) { // S2 cell
					b_u_32(index)   =  U(cell->neighbor(0)->neighbor_index(3)) - u0;  
					index++; 
                }

                if (cell->neighbor(3)->neighbor_index(1) != -1) { // S4 cell
					b_u_32(index)   =  U(cell->neighbor(3)->neighbor_index(1)) - u0;
				}

				CLS_R32[c].solve(b_u_32, d_u_32, u_coeff_32);
			}
            
            // If stencil is not available fallback to third order centered stencil
            
            else {
				u_coeff_32 = u_coeff_3; 
            }  
            
            // =====================================================================
            // r = 3 stencil 3
            // =====================================================================
            
            if (is_admissible_R33[c]) {
            
                // Form constraint RHS

				d_u_33(0)  = U(cell->neighbor_index(0)) - u0;                // P1 neighbor 
				d_u_33(1)  = U(cell->neighbor_index(1)) - u0;                // P3 neighbor
				d_u_33(2)  = U(cell->neighbor_index(2)) - u0;                // P4 neighbor
				d_u_33(3)  = U(cell->neighbor(2)->neighbor_index(2)) - u0;   // S7 neighbor
				
                if (cell->neighbor(1)->neighbor_index(2) != -1  &&
                    cell->neighbor(2)->neighbor_index(0) != -1) {
					ROWS = 2; 
                }

                else {
					ROWS = 1; 
                }

				b_u_33.reinit(ROWS);

				index = 0; 

                if (cell->neighbor(1)->neighbor_index(2) != -1) { // S2 cell
					b_u_33(index)   =  U(cell->neighbor(1)->neighbor_index(2)) - u0;  
					index++; 
                }

                if (cell->neighbor(2)->neighbor_index(0) != -1) { // S4 cell
					b_u_33(index)   =  U(cell->neighbor(2)->neighbor_index(0)) - u0;
                }

				CLS_R33[c].solve(b_u_33, d_u_33, u_coeff_33); 
            }
            
            // If stencil is not available fallback to third order centered stencil
            
            else {
				u_coeff_33 = u_coeff_3; 
            }
            
			// =====================================================================
            // r = 3 stencil 4
            // =====================================================================
            
            if (is_admissible_R34[c]) {
            
                // Form constraint RHS

				d_u_34(0)  = U(cell->neighbor_index(3)) - u0;                // P2 neighbor 
				d_u_34(1)  = U(cell->neighbor_index(1)) - u0;                // P3 neighbor
				d_u_34(2)  = U(cell->neighbor_index(2)) - u0;                // P4 neighbor
				d_u_34(3)  = U(cell->neighbor(1)->neighbor_index(1)) - u0;   // S5 neighbor			

                if (cell->neighbor(3)->neighbor_index(1) != -1  &&
                    cell->neighbor(2)->neighbor_index(1) != -1) {
					ROWS = 2; 
                }

                else {
					ROWS = 1; 
                }

				b_u_34.reinit(ROWS);

				index = 0; 

                if (cell->neighbor(3)->neighbor_index(1) != -1) { // S4 cell
					b_u_34(index)   =  U(cell->neighbor(3)->neighbor_index(1)) - u0;  
					index++; 
                }

                if (cell->neighbor(2)->neighbor_index(1) != -1) { // S6 cell
					b_u_34(index)   =  U(cell->neighbor(2)->neighbor_index(1)) - u0;
                }

				CLS_R34[c].solve(b_u_34, d_u_34, u_coeff_34);
            }
            
            // If stencil is not available fallback to third order centered stencil
            
            else {
				u_coeff_34 = u_coeff_3; 
            } 
            
			// =====================================================================
			// r = 2 stencil 1
			// ===================================================================== 

			b_u2(0)  = (U(cell->neighbor_index(0)) - u0);   // P1 neighbor 
			b_u2(1)  = (U(cell->neighbor_index(3)) - u0);   // P2 neighbor
			LU_R21[c].solve(b_u2, u_coeff_21);
            
            // =====================================================================
            // r = 2 stencil 2 
            // ===================================================================== 

			b_u2(0) = (U(cell->neighbor_index(3)) - u0);   // P2 neighbor
			b_u2(1) = (U(cell->neighbor_index(1)) - u0);   // P3 neighbor
			LU_R22[c].solve(b_u2, u_coeff_22);
		
			// =====================================================================
			// r = 2 stencil 3
			// =====================================================================

			b_u2(0) = (U(cell->neighbor_index(1)) - u0);   // P3 neighbor
			b_u2(1) = (U(cell->neighbor_index(2)) - u0);   // P4 neighbor                                           // P4 neighbor
			LU_R23[c].solve(b_u2, u_coeff_23); 
			
			// =====================================================================
			// r = 2 stencil 4
			// =====================================================================

			b_u2(0) = (U(cell->neighbor_index(2)) - u0);   // P4 neighbor
			b_u2(1) = (U(cell->neighbor_index(0)) - u0);   // P1 neighbor
			LU_R24[c].solve(b_u2, u_coeff_24);

        // Find the smoothness indicators 
        
		// =====================================================================
		// r = 4 Stencil 
		// =====================================================================
        
		IS_u(0)   = compute_fourth_order_smoothness_indicator(u_coeff_4, IS_constants[c], h);

		// =====================================================================
		// r = 3
		// =====================================================================
        
		IS_u(1)   = compute_third_order_smoothness_indicator(u_coeff_3, IS_constants[c], h); 
        
		// =====================================================================
		// r = 3  stencil 1
		// =====================================================================
        
		IS_u(2)   = compute_third_order_smoothness_indicator(u_coeff_31, IS_constants[c], h); 
        
        
		// =====================================================================
		// r = 3  stencil 2
		// =====================================================================
        
		IS_u(3)   = compute_third_order_smoothness_indicator(u_coeff_32, IS_constants[c], h); 
        
		// =====================================================================
		// r = 3  stencil 3
		// =====================================================================
        
		IS_u(4)   = compute_third_order_smoothness_indicator(u_coeff_33, IS_constants[c], h); 
        
		// =====================================================================
		// r = 3  stencil 4
		// =====================================================================
        
		IS_u(5)   = compute_third_order_smoothness_indicator(u_coeff_34, IS_constants[c], h); 
        
		// =====================================================================
		// r = 2  stencil 1
		// =====================================================================
		
		IS_u(6)   = compute_second_order_smoothness_indicator(u_coeff_21, IS_constants[c], h); 
		
		// =====================================================================
		// r = 2  stencil 2
		// =====================================================================
		
		IS_u(7)   = compute_second_order_smoothness_indicator(u_coeff_22, IS_constants[c], h); 
		
		// =====================================================================
		// r = 2  stencil 3
		// =====================================================================
		
		IS_u(8)   = compute_second_order_smoothness_indicator(u_coeff_23, IS_constants[c], h); 
		
		// =====================================================================
		// r = 2  stencil 4
		// =====================================================================
		
		IS_u(9)   = compute_second_order_smoothness_indicator(u_coeff_24, IS_constants[c], h); 
		
		// Combine the polynomials 
		
		sum_u = 0.0; sum_gamma = 0.0; 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			
			w_u(j)   = gamma(j)/(std::pow((IS_u(j) + epsilon), p));
			
			sum_u += w_u(j);
			
			sum_gamma += gamma(j); 
		}
		
		// Normalize the weights 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			w_u(j) = w_u(j)/sum_u; 
			gamma(j) = gamma(j)/sum_gamma; 
		}
		
		// Density 
		
		coeffs_U[c](0) = u0; 
		
		coeffs_U[c](1) = (w_u(0)/gamma(0))*(u_coeff_4(0) - 
												
												gamma(1)*u_coeff_3(0)  - 
												gamma(2)*u_coeff_31(0) - 
												gamma(3)*u_coeff_32(0) - 
												gamma(4)*u_coeff_33(0) - 
												gamma(5)*u_coeff_34(0) - 
												gamma(6)*u_coeff_21(0) -
												gamma(7)*u_coeff_22(0) -
												gamma(8)*u_coeff_23(0) -
												gamma(9)*u_coeff_24(0) ) + 
												
												w_u(1)*u_coeff_3(0) +
												w_u(2)*u_coeff_31(0) + 
												w_u(3)*u_coeff_32(0) +
												w_u(4)*u_coeff_33(0) +
												w_u(5)*u_coeff_34(0) +
												w_u(6)*u_coeff_21(0) +
												w_u(7)*u_coeff_22(0) +
												w_u(8)*u_coeff_23(0) +
												w_u(9)*u_coeff_24(0) ;// u_x
		
		coeffs_U[c](2) = (w_u(0)/gamma(0))*(u_coeff_4(1) - 
												
												gamma(1)*u_coeff_3(1) - 
												gamma(2)*u_coeff_31(1) - 
												gamma(3)*u_coeff_32(1) -
												gamma(4)*u_coeff_33(1) -
												gamma(5)*u_coeff_34(1) - 
												gamma(6)*u_coeff_21(1) -
												gamma(7)*u_coeff_22(1) -
												gamma(8)*u_coeff_23(1) -
												gamma(9)*u_coeff_24(1) ) + 
												
												w_u(1)*u_coeff_3(1) +
												w_u(2)*u_coeff_31(1) +
												w_u(3)*u_coeff_32(1) +
												w_u(4)*u_coeff_33(1) +
												w_u(5)*u_coeff_34(1) +
												w_u(6)*u_coeff_21(1) +
												w_u(7)*u_coeff_22(1) +
												w_u(8)*u_coeff_23(1) +
												w_u(9)*u_coeff_24(1) ; // u_x

		coeffs_U[c](3) = (w_u(0)/gamma(0))*(u_coeff_4(2) - 
												
												gamma(1)*u_coeff_3(2) - 
												gamma(2)*u_coeff_31(2) - 
												gamma(3)*u_coeff_32(2) - 
												gamma(4)*u_coeff_33(2) - 
												gamma(5)*u_coeff_34(2) ) + 
												
												w_u(1)*u_coeff_3(2) +
												w_u(2)*u_coeff_31(2) +
												w_u(3)*u_coeff_32(2) + 
												w_u(4)*u_coeff_33(2) + 
												w_u(5)*u_coeff_34(2) ; // u_xx
												
		coeffs_U[c](4) = (w_u(0)/gamma(0))*(u_coeff_4(3) - 
												
												gamma(1)*u_coeff_3(3) - 
												gamma(2)*u_coeff_31(3) - 
												gamma(3)*u_coeff_32(3) - 
												gamma(4)*u_coeff_33(3) - 
												gamma(5)*u_coeff_34(3) ) + 
												
												w_u(1)*u_coeff_3(3) +
												w_u(2)*u_coeff_31(3) + 
												w_u(3)*u_coeff_32(3) + 
												w_u(4)*u_coeff_33(3) + 
												w_u(5)*u_coeff_34(3) ; // u_yy
		
		coeffs_U[c](5) = (w_u(0)/gamma(0))*(u_coeff_4(4) - 
												
												gamma(1)*u_coeff_3(4) - 
												gamma(2)*u_coeff_31(4) - 
												gamma(3)*u_coeff_32(4) - 
												gamma(4)*u_coeff_33(4) - 
												gamma(5)*u_coeff_34(4) ) + 
												
												w_u(1)*u_coeff_3(4) +
												w_u(2)*u_coeff_31(4) + 
												w_u(3)*u_coeff_32(4) + 
												w_u(4)*u_coeff_33(4) + 
												w_u(5)*u_coeff_34(4) ; // u_xy
		
		coeffs_U[c](6) = (w_u(0)/gamma(0))*(u_coeff_4(5)); // u_xxx

		coeffs_U[c](7) = (w_u(0)/gamma(0))*(u_coeff_4(6)); // u_yyy

		coeffs_U[c](8) = (w_u(0)/gamma(0))*(u_coeff_4(7)); // u_xxy

		coeffs_U[c](9) = (w_u(0)/gamma(0))*(u_coeff_4(8)); // u_xyy

            
		} // End of center block loop 
        
    } // End of cell loop 
    
} // End of function 

