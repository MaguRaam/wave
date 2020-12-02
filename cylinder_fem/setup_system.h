 template <int dim>
  void WaveEquation<dim>::setup_system ()
  {
  		double x_min = -1.0,  
    			 x_max = 1.0,  
    			 y_min = -1.0,  
				 y_max = 1.0;
				 
  		Point<dim,double> min(x_min,y_min),
								max(x_max,y_max);
		std::vector<unsigned int> numberofcells{200,200};
  
  		GridGenerator::subdivided_hyper_rectangle (triangulation, numberofcells, min, max);
       
		
		//graphical output of the generated grid:
		
		std::ofstream out("grid.eps");
  		GridOut       grid_out;
      grid_out.write_eps(triangulation, out);
      std::cout << "Grid written to grid-2.eps" << std::endl;
		std::cout << "Mesh info:" << std::endl
            << " dimension: " << 2 << std::endl
            << " no. of cells: " << triangulation.n_active_cells() << std::endl;
		std::cout << "Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl;

    dof_handler.distribute_dofs (fe);

    std::cout << "Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl
              << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
     
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    sparsity_pattern.copy_from (dsp);

	 
    mass_matrix.reinit (sparsity_pattern);
    laplace_matrix.reinit (sparsity_pattern);
    matrix_u.reinit (sparsity_pattern);
    matrix_v.reinit (sparsity_pattern);
 
	 QGauss<dim> quadrature_formula(fe.degree + 1);
		FEValues<dim> fe_values(fe,
				  quadrature_formula,
				  update_values | 
				  update_gradients |
				  update_quadrature_points| 
				  update_JxW_values);
				  
		const unsigned int 			dofs_per_cell = fe.dofs_per_cell;
		unsigned int       			num_quad_pts = quadrature_formula.size();
		FullMatrix<double> 			Mlocal (dofs_per_cell, dofs_per_cell);
		FullMatrix<double> 			Llocal (dofs_per_cell, dofs_per_cell);
		std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

for (const auto &cell : dof_handler.active_cell_iterators())
  {
    Llocal = 0.;
    Mlocal = 0.;
    
    fe_values.reinit(cell);
    for (unsigned int q_index = 0; q_index < num_quad_pts; ++q_index)
      {
        const double current_a =
          a<dim>(fe_values.quadrature_point(q_index));  
          
        const double current_rho =
          rho<dim>(fe_values.quadrature_point(q_index));    
          
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {Llocal(i, j) +=
                (current_a *              // a(x_q)
                 fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx
                 
               Mlocal(i,j)+=current_rho*fe_values.shape_value(i, q_index)*fe_values.shape_value(j, q_index)*fe_values.JxW(q_index);  
              }
          }
      }
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {laplace_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            Llocal(i, j));
          mass_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            Mlocal(i, j));}
         
      }
  }
 
 	                               
                                          
                                          
                                          
                                          
                                          
 	 wave_speed.reinit (dof_handler.n_dofs());
    solution_u.reinit (dof_handler.n_dofs());
    solution_v.reinit (dof_handler.n_dofs());
    old_solution_u.reinit (dof_handler.n_dofs());
    old_solution_v.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    constraints.close ();
  }

