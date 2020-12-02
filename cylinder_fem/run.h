template <int dim>
  void WaveEquation<dim>::run ()
  {
    setup_system();

    VectorTools::project (dof_handler, constraints, QGauss<dim>(fe.degree + 1),
                          InitialValuesU<dim>(),
                          old_solution_u);
    VectorTools::project (dof_handler, constraints, QGauss<dim>(fe.degree + 1),
                          InitialValuesV<dim>(),
                          old_solution_v);
    VectorTools::project (dof_handler, constraints, QGauss<dim>(fe.degree + 4),
                          InitialValueswave_speed<dim>(),
                          wave_speed);                      
                          

    Vector<double> tmp (solution_u.size());
    Vector<double> forcing_terms (solution_u.size());

	  unsigned int c = 0;
		  
      
		  
		  
    for (; time<=5; time+=time_step, ++timestep_number)
      {
        std::cout << "Time step " << timestep_number
                  << " at t=" << time
                  << std::endl;

        mass_matrix.vmult (system_rhs, old_solution_u);

        mass_matrix.vmult (tmp, old_solution_v);
        system_rhs.add (time_step, tmp);

        laplace_matrix.vmult (tmp, old_solution_u);
        system_rhs.add (-theta * (1-theta) * time_step * time_step, tmp);

        RightHandSide<dim> rhs_function;
        rhs_function.set_time (time);
        VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(fe.degree + 1),
                                             rhs_function, tmp);
        forcing_terms = tmp;
        forcing_terms *= theta * time_step;

        rhs_function.set_time (time-time_step);
        VectorTools::create_right_hand_side (dof_handler, QGauss<dim>(fe.degree + 1),
                                             rhs_function, tmp);

        forcing_terms.add ((1-theta) * time_step, tmp);

        system_rhs.add (theta * time_step, forcing_terms);

        {
          BoundaryValuesU<dim> boundary_values_u_function;
          boundary_values_u_function.set_time (time);

          std::map<types::global_dof_index,double> boundary_values;
          VectorTools::interpolate_boundary_values (dof_handler,
                                                    0,
                                                    boundary_values_u_function,
                                                    boundary_values);
          

          matrix_u.copy_from (mass_matrix);
          matrix_u.add (theta * theta * time_step * time_step, laplace_matrix);
          MatrixTools::apply_boundary_values (boundary_values,
                                              matrix_u,
                                              solution_u,
                                              system_rhs);
        }
        solve_u ();


        laplace_matrix.vmult (system_rhs, solution_u);
        system_rhs *= -theta * time_step;

        mass_matrix.vmult (tmp, old_solution_v);
        system_rhs += tmp;

        laplace_matrix.vmult (tmp, old_solution_u);
        system_rhs.add (-time_step * (1-theta), tmp);

        system_rhs += forcing_terms;

        {
          BoundaryValuesV<dim> boundary_values_v_function;
          boundary_values_v_function.set_time (time);

          std::map<types::global_dof_index,double> boundary_values;
          VectorTools::interpolate_boundary_values (dof_handler,
                                                    0,
                                                    boundary_values_v_function,
                                                    boundary_values);
                                                    
         
          matrix_v.copy_from (mass_matrix);
          MatrixTools::apply_boundary_values (boundary_values,
                                              matrix_v,
                                              solution_v,
                                              system_rhs);
        }
        solve_v ();
		 
        if(c%100 == 0) output_results ();
        c++;

        std::cout << "   Total energy: "
                  << (mass_matrix.matrix_norm_square (solution_v) +
                      laplace_matrix.matrix_norm_square (solution_u)) / 2
                  << std::endl;

        old_solution_u = solution_u;
        old_solution_v = solution_v;
      }
  }
 
