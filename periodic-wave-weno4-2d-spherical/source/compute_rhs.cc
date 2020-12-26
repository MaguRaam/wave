#include "../include/Weno432.h"


// Compute the rhs vectors

void Weno4_2D::compute_rhs() 
{


    copy_data();
    reconstruct();
    
    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;

    auto cell = dof_handler.begin_active();
    auto endc = dof_handler.end();

    Point<2> face_quadrature_point_1;     
    Point<2> periodic_face_quadrature_point_1;
    Point<2> face_quadrature_point_2;     
    Point<2> periodic_face_quadrature_point_2;

    unsigned int neighbor_index;

    double V_c;

    double nx1, ny1;    
    double nx2, ny2;
    double j_w1, j_w2;

    double GL1_x, GR1_x;
    double GL1_y, GR1_y;

    double GL2_x, GR2_x;
    double GL2_y, GR2_y;

    double F;

    unsigned int global_face_index;   
    
    unsigned int n_faces = triangulation.n_active_faces(); 
 
    std::vector<double> Flux1(n_faces);
    std::vector<double> Flux2(n_faces);
    
    std::vector< bool > did_not_compute_flux_for_the_face(n_faces); 

    for (unsigned int f = 0; f < n_faces; f++) {
        did_not_compute_flux_for_the_face[f] = true; 
    }



    for (unsigned int c = 0; cell != endc; ++cell, ++c) 
    {
        
        if (c < no_cells_per_block) 
        {

            V_c = Cell[c].measure();

            rhs1(c) = V(c);
            rhs2(c) = 0.0;
        
            for (unsigned int f = 0; f < faces_per_cell; ++f) 
            {

                global_face_index = cell->face_index(f);

                j_w1 = Cell[c].face_jxw_1(f);
                j_w2 = Cell[c].face_jxw_2(f);

                if(did_not_compute_flux_for_the_face[global_face_index])
                {

                    double length = 4.0; //periodic length

                    nx1 = Cell[c].nx1(f); ny1 = Cell[c].ny1(f);
                    nx2 = Cell[c].nx2(f); ny2 = Cell[c].ny2(f);

                    face_quadrature_point_1 = Cell[c].face_quadrature_point1(f); 
                    face_quadrature_point_2 = Cell[c].face_quadrature_point2(f); 

                    //Left face:
                    GL1_x = evaluate_weno_gradient_x(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_1);
                    GL2_x = evaluate_weno_gradient_x(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_2);

                    GL1_y = evaluate_weno_gradient_y(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_1);
                    GL2_y = evaluate_weno_gradient_y(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_2);
                    
                    // Get the right state values
                    neighbor_index = cell->neighbor_index(f);

                    if (neighbor_index >= 2*no_cells_per_block && neighbor_index < 3*no_cells_per_block) {
                        
                        // Bottom 
                        
                        periodic_face_quadrature_point_1 = face_quadrature_point_1;           
                        periodic_face_quadrature_point_2 = face_quadrature_point_2;                   
                        
                        periodic_face_quadrature_point_1(1) =  periodic_face_quadrature_point_1(1) + length;
                        periodic_face_quadrature_point_2(1) =  periodic_face_quadrature_point_2(1) + length;
                        
                        GR1_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 2*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 2*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 2*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 2*no_cells_per_block], periodic_face_quadrature_point_2);

                        GR1_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 2*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 2*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 2*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 2*no_cells_per_block], periodic_face_quadrature_point_2);
                        
                    }

                    else if (neighbor_index >= 4*no_cells_per_block && neighbor_index < 5*no_cells_per_block) {

                        // Left 
                        
                        periodic_face_quadrature_point_1 = face_quadrature_point_1;           
                        periodic_face_quadrature_point_2 = face_quadrature_point_2;                   
                        
                        periodic_face_quadrature_point_1(0) =  periodic_face_quadrature_point_1(0) + length;
                        periodic_face_quadrature_point_2(0) =  periodic_face_quadrature_point_2(0) + length;
                        
                        GR1_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 4*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 4*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 4*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 4*no_cells_per_block], periodic_face_quadrature_point_2);
                        
                        GR1_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 4*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 4*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 4*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 4*no_cells_per_block], periodic_face_quadrature_point_2);
                    }

                    else if (neighbor_index >= 5*no_cells_per_block && neighbor_index < 6*no_cells_per_block) {

                        // Right 
                        
                        periodic_face_quadrature_point_1 = face_quadrature_point_1;           
                        periodic_face_quadrature_point_2 = face_quadrature_point_2;                   
                        
                        periodic_face_quadrature_point_1(0) =  periodic_face_quadrature_point_1(0) - length;
                        periodic_face_quadrature_point_2(0) =  periodic_face_quadrature_point_2(0) - length;
                        
                        GR1_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 5*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 5*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 5*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 5*no_cells_per_block], periodic_face_quadrature_point_2);
                        
                        GR1_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 5*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 5*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 5*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 5*no_cells_per_block], periodic_face_quadrature_point_2);
                        
                    }

                    else if (neighbor_index >= 7*no_cells_per_block && neighbor_index < 8*no_cells_per_block) {

                        // Top 
                        
                        periodic_face_quadrature_point_1 = face_quadrature_point_1;           
                        periodic_face_quadrature_point_2 = face_quadrature_point_2;                   
                        
                        periodic_face_quadrature_point_1(1) =  periodic_face_quadrature_point_1(1) - length;
                        periodic_face_quadrature_point_2(1) =  periodic_face_quadrature_point_2(1) - length;
                        
                        GR1_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 7*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 7*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_x = evaluate_weno_gradient_x(coeffs_U[neighbor_index - 7*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 7*no_cells_per_block], periodic_face_quadrature_point_2);
                        
                        GR1_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 7*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 7*no_cells_per_block], periodic_face_quadrature_point_1);
                        
                        GR2_y = evaluate_weno_gradient_y(coeffs_U[neighbor_index - 7*no_cells_per_block], 
                                                         WENO_poly_consts[neighbor_index - 7*no_cells_per_block], periodic_face_quadrature_point_2);
                          
                    }

                    else {
                    
                        GR1_x = evaluate_weno_gradient_x(coeffs_U[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1);
                        
                        GR2_x = evaluate_weno_gradient_x(coeffs_U[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2);

                        GR1_y = evaluate_weno_gradient_y(coeffs_U[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1);
                        
                        GR2_y = evaluate_weno_gradient_y(coeffs_U[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2);
                    }

                    Flux1[global_face_index] = gradient_flux(GL1_x,GR1_x,GL1_y,GR1_y,nx1,ny1);
                    Flux2[global_face_index] = gradient_flux(GL2_x,GR2_x,GL2_y,GR2_y,nx2,ny2);
                                        
                    did_not_compute_flux_for_the_face[global_face_index] = false;
                }

                else 
                {
                    
                    Flux1[global_face_index] *= -1.0; 
                    Flux2[global_face_index] *= -1.0;
                }

                F = j_w1 * Flux1[global_face_index] + j_w2 * Flux2[global_face_index];  

                // Add it to the rhs vectors

                rhs2(c) += (1.0/V_c)*F;


            }

            rhs2(c) = C*C*rhs2(c);


        } // End of center block loop

    }   // End of cell loop 
         
}
