template <int dim>
  void WaveEquation<dim>::solve_u ()
  {
    SolverControl           solver_control (1000, 1e-8*system_rhs.l2_norm());
    SolverCG<>              cg (solver_control);

    cg.solve (matrix_u, solution_u, system_rhs,
              PreconditionIdentity());

    std::cout << "   u-equation: " << solver_control.last_step()
              << " CG iterations."
              << std::endl;
  }


  template <int dim>
  void WaveEquation<dim>::solve_v ()
  {
    SolverControl           solver_control (1000, 1e-8*system_rhs.l2_norm());
    SolverCG<>              cg (solver_control);

    cg.solve (matrix_v, solution_v, system_rhs,
              PreconditionIdentity());

    std::cout << "   v-equation: " << solver_control.last_step()
              << " CG iterations."
              << std::endl;
  }
