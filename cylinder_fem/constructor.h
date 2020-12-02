template <int dim>
  WaveEquation<dim>::WaveEquation () :
    fe (1),
    dof_handler (triangulation),
    time_step (0.001),
    time (time_step),
    timestep_number (1),
    theta (0.5)
  {}

//See .geo file:
//Length of inner circle=1.5707
//No of grid=100
