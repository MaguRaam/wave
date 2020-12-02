#include "fv_tools.h"
#include "ode.h"

/*Solve for u(x,t),x ∈ Ω([−1,1]^2),t ∈[0,T]
  ∂^2u/∂^2t−∇^2u= 0   x∈Ω

  Boundary condition:
  u(x,t) = 0      x∈∂Ω

  Initial condition:
  u(x,0) =sin(πkx)sin(πly)  x∈Ω

  ∂u/∂t(x,0) = 0      x∈Ω
*/



template <int dim> class standing_wave {
public:
  standing_wave(int refine);
  Vector<double> initialize();
  void operator()(const Vector<double> &x, Vector<double> &dxdt, double /*t*/);
  void write_data(int timestep_number);

private:
  Triangulation<dim> triangulation;
  DoFHandler<dim> dof_handler;
  Reconstruct<dim> reconstruct;
  CellLaplace<dim> laplace;
  CellAverage<dim, 4> cell_average;

  unsigned int n_cells;

  Vector<double> u, u_exact,u_error;
  GlobalPolynomial<dim> u_poly;
  BoundaryCondition<dim> u_bcs;
};

// constructor:
template <int dim>
standing_wave<dim>::standing_wave(int refine) : dof_handler(triangulation) {
  GridGenerator::hyper_cube(triangulation, -1., 1., false);
  triangulation.refine_global(refine);
  //dealii::GridTools::distort_random(0.1, triangulation, true);

  // write mesh file:
  std::ofstream out("plot/grid.vtk");
  GridOut grid_out;
  grid_out.write_vtk(triangulation, out);

  // enumerate finite volume cells:
  dof_handler.distribute_dofs(FE_DGQ<dim>{0});
  n_cells = dof_handler.n_dofs();

  std::cout << "n_cells = " << n_cells << "\n";

  // create polynomial matrix for each cell:
  reconstruct.create_matrix(dof_handler);

  // resize solution and polynomial:
  u.reinit(n_cells);
  u_exact.reinit(n_cells);
  u_poly.resize(n_cells);
  u_error.reinit(n_cells);
  

  // set dirchlet bcs u = 0:
  BoundaryCondition<dim> u_bcs;  
  for (const auto &cell : dof_handler.active_cell_iterators()){
    if (cell->at_boundary()){
      for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary()) u_bcs[cell->face(f)] = 0.0;
    }
  }

}

// initialize state x:
template <int dim> Vector<double> standing_wave<dim>::initialize() {
  Vector<double> x(2 * n_cells);
  x = 0.0;
  int k(2);
  auto uf = [k](auto p) { return sin(M_PI * k * p(0)) * sin(M_PI * k * p(1)) * sin(M_PI * k * p(2)); };

  unsigned int i(0);
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    x[i] = cell_average(uf, cell);
    i++;
  }
  return x;
}

// write data:
template <int dim> void standing_wave<dim>::write_data(int timestep_number) {
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(u, "u", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(u_exact, "u_exact", DataOut<dim>::type_dof_data);
  data_out.add_data_vector(u_error, "u_error", DataOut<dim>::type_dof_data);
  data_out.build_patches();

  const std::string filename = "plot/solution-" +
                               Utilities::int_to_string(timestep_number, 3) +
                               ".vtk";
  DataOutBase::VtkFlags vtk_flags;
  vtk_flags.compression_level =
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed;
  data_out.set_flags(vtk_flags);
  std::ofstream output(filename);
  data_out.write_vtk(output);
}

// dx/dt = f(x) or du/dt = v , dv/dt = L(u)
template <int dim>
void standing_wave<dim>::operator()(const Vector<double> &x,
                                    Vector<double> &dxdt, double t) {

  // du/dt = v:
  using std::copy;
  copy(x.begin() + n_cells, x.end(), dxdt.begin());

  // extract u:
  copy(x.begin(), x.begin() + n_cells, u.begin());

  // Reconstruct polynomial for u:
  reconstruct(dof_handler, u, u_bcs, u_poly);

  // dv/dt = L(u):
  unsigned int i = n_cells;
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    dxdt[i] = laplace(cell, u_poly);
    i++;
  }

  // exact solution:
  int k(2);
  double omega(3 * M_PI * M_PI * k * k);
  auto uf_exact = [k, t, omega](auto p) {
    return (cos(sqrt(omega) * t)) * sin(M_PI * k * p(0)) * sin(M_PI * k * p(1)) * sin(M_PI * k * p(2));
  };

  i = 0;
  for (const auto &cell : dof_handler.active_cell_iterators()) {
    u_exact[i] = cell_average(uf_exact, cell);
    u_error[i] = abs(u_exact[i] - u[i]);
    i++;
  }

  std::cout.flags( std::ios::dec | std::ios::scientific );
  std::cout.precision(5);
  std::cout<<"t = "<<t<<"  Linfty error = "<<linfty_error<dim>(u,uf_exact,dof_handler)<<"\n";
  std::cout<<"t = "<<t<<"  L2 error = "<<l2_error<dim>(u,uf_exact,dof_handler)<<"\n";
}

using RK4 = runge_kutta4<Vector<double>>;

void convergence(int refinement){
  standing_wave<3> wave(refinement);

  RK4 rk4;
  const int steps = 1000;
  const double dt = 0.001;

  // initialize:
  Vector<double> x(wave.initialize());
  for (size_t n = 0; n < steps; ++n)
  {
        rk4.do_step(wave, x, n * dt, dt);
        std::cout << n * dt << '\n';
        if (n%100 == 0) wave.write_data(n);
  }

}



int main() {
    convergence(4); //8*8
  return 0;
}


 